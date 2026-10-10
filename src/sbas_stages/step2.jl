"""
    step2(
        data::DataInfo,
        chunks::Vector{ChunkInfo},
        tile::TileInfo,
        sys::SystemInfo,
        Binv::Matrix{Float32},
        wavelength::Float64,
        tkp1mtk::Vector{Float32}
    )

Runs the inversion at all pixels to derive the deformation time series.

# Notes
`Binv` is the pseudo-inverse of the design matrix. It maps the stack of 
interferometric phases [rad] onto the velocity of each epoch interval 
[rad/day], which is then integrated with `tkp1mtk` and scaled to [cm].
"""
function step2(
    data::DataInfo,
    chunks::Vector{ChunkInfo},
    tile::TileInfo,
    sys::SystemInfo,
    Binv::Matrix{Float32},
    wavelength::Float64,
    tkp1mtk::Vector{Float32}
)

    # Conversion from phase to cm
    coef_def = Float32(wavelength / (4 * pi))

    # Buffers
    chunk_buffer = Matrix{Float32}(undef, chunks[1].numel, data.nint)
    dhist_buffer = Matrix{Float32}(undef, chunks[1].numel, data.nacq - 1)

    # Open time series file 
    wio = open(tile.path, "a")
    try
        for (i, chunk_info) in enumerate(chunks)
            println("Processing chunk $i of $(length(chunks))")

            # Read chunk into buffer
            chunk = view(chunk_buffer, 1:chunk_info.numel, :)
            open(chunk_info.path, "r") do rio
                read!(rio, chunk)
            end

            # Parallel processing
            pix_per_thread = ceil(Int, chunk_info.numel / sys.threads)
            @threads for tid in 1:sys.threads

                # Thread slice 
                row1 = pix_per_thread * (tid - 1) + 1
                row2 = min(pix_per_thread * tid, chunk_info.numel)

                # Only proceed if needed 
                # <<< in case there are less pixels than threads >>>
                if row1 <= row2

                    # Inversion
                    for row in row1:row2

                        # Epochs velocity
                        d = view(chunk, row, :)
                        v = view(dhist_buffer, row, :)
                        mul!(v, Binv, d)

                        # Integrate
                        v[1] *= tkp1mtk[1]
                        for k in 2:data.nacq-1
                            v[k] = v[k-1] + v[k] * tkp1mtk[k]
                        end
                    end
                end
            end

            # Write out deformation time series
            dhist = view(dhist_buffer, 1:chunk_info.numel, :)
            dhist .*= coef_def  # Convert to deformation [cm]
            for t in 1:(data.nacq-1)
                element_offset = (t - 1) * data.npix + (chunk_info.p1 - 1)
                byte_offset = element_offset * sizeof(Float32)
                seek(wio, byte_offset)
                write(wio, view(dhist, :, t))
            end
        end
    finally
        close(wio)
    end
end

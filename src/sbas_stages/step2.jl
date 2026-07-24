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
    coef = wavelength / (4 * pi)

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
            @threads for k in 1:sys.threads

                # Thread slice 
                row1 = pix_per_thread * (k - 1) + 1
                row2 = min(pix_per_thread * k, chunk_info.numel)

                # Only proceed if needed (in case there are less pixels than threads)
                if row1 <= row2

                    # Inversion (solve for velocity then integrate)
                    for row in row1:row2
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

            # Write out deformation
            vhist = view(dhist_buffer, 1:chunk_info.numel, :)
            vhist .*= coef  # Convert to deformation in cm
            for t in 1:(data.nacq-1)
                element_offset = (t - 1) * data.npix + (chunk_info.p1 - 1)
                byte_offset = element_offset * sizeof(Float32)
                seek(wio, byte_offset)
                write(wio, view(vhist, :, t))
            end
        end
    finally
        close(wio)
    end

end

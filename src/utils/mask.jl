"""
    mask_unw_igrams(
        wid::Int,
        len::Int,
        mask_path::String,
        intlist::String
    )

Apply a mask to a stack of unwrapped interferograms.

# Positional arguments
- `wid`: `Int` width of unwrapped interferograms (in pixels)
- `len`: `Int` length of unwrapped interferograms (in pixels)
- `mask_path`: `String` path to mask
- `intlist`: `String` path to file containing the list of interferograms

# Returns
Masked interferograms in the `masked` directory
"""
function mask_unw_igrams(
    wid::Int,
    len::Int,
    mask_path::String,
    intlist::String
)
    # Pre-allocate buffers
    raw = Matrix{Float32}(undef, 2 * wid, len)
    amp = Matrix{Float32}(undef, len, wid)
    psi = Matrix{Float32}(undef, len, wid)
    mask_float = Matrix{Float32}(undef, len, wid)

    # Read the mask
    read_unw_igram!(raw, amp, mask_float, mask_path)
    valid = mask_float .> Float32(0.5)

    # Read list of files
    file_paths = Vector{String}()
    open(intlist, "r") do io
        for line in eachline(io)
            if !isempty(strip(line))
                # Check the file exists
                if isfile(line)
                    push!(file_paths, line)
                else
                    error("Could not find: $(line)")
                end
            end
        end
    end

    # Create the output directory
    dirin = dirname(file_paths[1])
    dirout = joinpath(dirin, "masked")
    if isdir(dirout)
        rm(dirout, recursive=true)
    end
    mkdir(dirout)

    # Mask interferograms
    for file_path in file_paths
        read_unw_igram!(raw, amp, psi, file_path)
        @inbounds @simd for k in eachindex(amp)
            if !valid[k]
                amp[k] = 0.0f0
                psi[k] = 0.0f0
            end
        end
        pathout = joinpath(dirout, basename(file_path))
        write_unw_igram!(raw, pathout, amp, psi)
    end
end

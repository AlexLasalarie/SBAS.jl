"""
    unwrapping_error(
        wid::Int,
        len::Int,
        intlist::String;
        name_out::String="unwrapping_error"
    )

Compute the total unwrapping error for each interferogram.

# Positional Arguments
- `wid`: `Int` width of unwrapped interferograms (in pixels)
- `len`: `Int` length of unwrapped interferograms (in pixels)
- `intlist`: `String` path to file containing list of interferograms

# Optional Arguments
- `name_out`: `String` name of output file (default is `unwrapping_error`)

# Returns
A file containing the total unwrapping error for each interferogram.
"""
function unwrapping_error(
    wid::Int,
    len::Int,
    intlist::String;
    name_out::String="unwrapping_error"
)
    # Pre-allocate buffers
    raw = Matrix{Float32}(undef, 2 * wid, len)
    amp = Matrix{Float32}(undef, len, wid)
    psi = Matrix{Float32}(undef, len, wid)

    # Read list of files
    paths = Vector{String}()
    open(intlist, "r") do io
        for line in eachline(io)
            if !isempty(strip(line))
                # Check the file exists
                if isfile(line)
                    push!(paths, line)
                else
                    error("Could not find: $(line)")
                end
            end
        end
    end
    nint = length(paths)
    dirpath = dirname(paths[1])

    # Compute total unwrapping error
    println("Computing unwrapping error...")
    errors = Vector{Float32}(undef, nint)
    for (k, path) in enumerate(paths)
        read_unw_igram!(raw, amp, psi, path)
        errors[k] = unw_err(psi)
    end

    # Output the error to file
    path_out = joinpath(dirpath, name_out)
    println("Outputing to: $path_out")
    open(path_out, "w") do io
        for (path, error) in zip(paths, errors)
            println(io, "$path $(error)")
        end
    end
    return path_out
end

function unw_err(psi::Matrix{Float32})
    len, wid = size(psi)
    unw_error = Float32(0.0)
    pi_f32 = Float32(pi)
    for col in 1:wid-1
        for row in 1:len-1
            a = psi[row, col]
            b = psi[row+1, col]
            c = psi[row+1, col+1]
            d = psi[row, col+1]
            for (p1, p2) in ((a, b), (b, c), (c, d), (d, a))
                dif = abs(p2 - p1)
                if dif > pi_f32
                    unw_error += dif
                end
            end
        end
    end
    return unw_error
end

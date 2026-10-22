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
    valid = BitMatrix(undef, len, wid)

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
        valid .= amp .> 0
        errors[k] = unw_err_edge(psi, valid)
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

"""
Edge-based unwrapping error.
"""
function unw_err_edge(
    psi::Matrix{Float32},
    isvalid::BitMatrix
)
    len, wid = size(psi)
    unwerr = zero(Float32)
    pi32 = Float32(pi)
    for j in 1:wid
        for i in 1:len
            val11 = psi[i, j]

            # Column difference
            if j + 1 <= wid     # skip if at right edge
                val12 = psi[i, j+1]
                if isvalid[i, j] && isvalid[i, j+1]
                    dif = abs(val12 - val11)
                    if dif > pi32
                        unwerr += dif
                    end
                end
            end

            # Row difference
            if i + 1 <= len     # skip if at bottom
                val21 = psi[i+1, j]
                if isvalid[i, j] && isvalid[i+1, j]
                    dif = abs(val21 - val11)
                    if dif > pi32
                        unwerr += dif
                    end
                end
            end
        end
    end
    return unwerr
end

"""
Cell-based unwrapping error, double counts interior points.
"""
function unw_err_cell(
    psi::Matrix{Float32},
    valid::BitMatrix
)
    len, wid = size(psi)
    unw_error = Float32(0.0)
    pi_f32 = Float32(pi)
    for col in 1:wid-1
        for row in 1:len-1
            # Values
            a = psi[row, col]
            b = psi[row+1, col]
            c = psi[row+1, col+1]
            d = psi[row, col+1]

            # Check validity
            va = valid[row, col]
            vb = valid[row+1, col]
            vc = valid[row+1, col+1]
            vd = valid[row, col+1]

            # Form (value, valid) pairs
            pairs = (
                ((a, va), (b, vb)),
                ((b, vb), (c, vc)),
                ((c, vc), (d, vd)),
                ((d, vd), (a, va))
            )

            # Update running sum
            for ((p1, v1), (p2, v2)) in pairs
                if v1 && v2
                    dif = abs(p2 - p1)
                    if dif > pi_f32
                        unw_error += dif
                    end
                end
            end
        end
    end
    return unw_error
end

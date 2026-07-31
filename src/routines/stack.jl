"""
    stack_igrams(
        wid::Int,
        len::Int,
        ref_col::Int,
        ref_row::Int,
        intlist::String;
        wavelength::Real=5.55
    )

Computes the linear deformation rate (in cm/year) using stacking as:
    v = lambda/(4*pi) * sum(psi)/sum(dt)

# Positional Arguments
- `wid`: `Int` width of unwrapped interferograms (in pixels)
- `len`: `Int` length of unwrapped interferograms (in pixels)
- `intlist`: `String` path to the file containing the list of files

# Optional Arguments
- `wavelength`: `Real` wavelength of sensor (in cm, default is `5.55`)

# Returns
`stack_rate.un`: file of standard `unwrapped interferogram` format containing 
the average amplitude (band 1) and linear deformation rate (in cm/year, band 2) 
"""
function stack_igrams(
    wid::Int,
    len::Int,
    ref_col::Int,
    ref_row::Int,
    intlist::String;
    wavelength::Real=5.55
)
    # Buffers
    raw = Matrix{Float32}(undef, 2 * wid, len)
    amp = Matrix{Float32}(undef, len, wid)
    psi = Matrix{Float32}(undef, len, wid)
    amp_sum = zeros(Float32, len, wid)
    psi_sum = zeros(Float32, len, wid)

    # Read list of files
    paths = parse_file(intlist, "s")[1]
    nint = length(paths)

    # Stack
    dtsum_int = Int(0)
    for path in paths
        println("Stacking $path")
        date1, date2 = extract_dates(path)
        dt = Int((date2 - date1).value)
        read_unw_igram!(raw, amp, psi, path)
        calibrate!(psi, ref_row, ref_col)
        if dt < 0
            psi .*= Float32(-1.0)
            dt = -dt
        end
        amp_sum .+= amp
        psi_sum .+= psi
        dtsum_int += dt
    end

    # Convert to deformation
    coef = Float32(wavelength) / (4 * pi)
    amp_sum ./= nint
    psi_sum .*= coef                        # in cm

    # Deformation rate (cm/year)
    dtsum = Float32(dtsum_int / 365.25)     # in years
    dirout = dirname(paths[1])
    path_out = joinpath(dirout, "stack_rate.un")
    write_unw_igram!(raw, path_out, amp_sum, psi_sum ./ dtsum)
end

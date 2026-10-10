"""
    stack_igrams(
        wid::Int,
        len::Int,
        ref_col::Int,
        ref_row::Int,
        intlist::String;
        wavelength::Real=5.55
    )

Computes the linear deformation rate [cm/year] using stacking as:
    v = lambda/(4*pi) * sum(psi)/sum(dt)

Also computes the linear deformation rate [cm/year] using linear SBAS as:
    v = lambda/(4*pi) * sum(psi*dt)/sum(dt^2)

Pairs named `date2_date1` (reversed chronological order) are negated before 
stacking, so either naming order may appear in `intlist`.

# Positional Arguments
- `wid`: `Int` width of unwrapped interferograms (in pixels)
- `len`: `Int` length of unwrapped interferograms (in pixels)
- `ref_col`: `Int` column location of the reference pixel
- `ref_row`: `Int` row location of the reference pixel
- `intlist`: `String` path to the file containing the list of files

# Optional Arguments
- `wavelength`: `Real` wavelength of sensor (in cm, default is `5.55`)

# Returns
`stack_rate.un`: file of standard `unwrapped interferogram` format containing 
the average amplitude (band 1) and the stacking deformation rate (in cm/year, 
band 2).

`sbas_rate.un`: file of standard `unwrapped interferogram` format containing 
the average amplitude (band 1) and the linear-SBAS deformation rate (in 
cm/year, band 2). This is the OLS fit through the origin over the 
interferogram pairs.
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
    psidt_sum = zeros(Float32, len, wid)

    # Read list of files
    paths = parse_file(intlist, "s")[1]
    nint = length(paths)

    # Stack
    dtsum_day = Float32(0.0)
    dt2sum_day = Float32(0.0)
    for path in paths
        println("Stacking $path")
        date1, date2 = extract_dates(path)
        dt = Int((date2 - date1).value)
        read_unw_igram!(raw, amp, psi, path)
        if dt < 0
            psi .*= Float32(-1.0)
            dt = -dt
        end
        calibrate!(psi, ref_row, ref_col)
        amp_sum .+= amp
        psi_sum .+= psi
        psidt_sum .+= psi .* dt
        dtsum_day += dt
        dt2sum_day += dt * dt
    end

    # Convert to deformation
    coef = Float32(wavelength) / (4 * pi)
    amp_sum ./= nint
    psi_sum .*= coef                            # in cm

    # Stacking deformation rate [cm/year]
    dtsum_year = Float32(dtsum_day / 365.25)     # in years
    dirout = dirname(paths[1])
    path_out = joinpath(dirout, "stack_rate.un")
    write_unw_igram!(raw, path_out, amp_sum, psi_sum ./ dtsum_year)

    # Linear SBAS deformation rate [cm/year]
    psidt_sum .*= coef              # in cm * day
    psidt_sum ./= dt2sum_day        # in (cm * day) / day^2 -> cm / day
    psidt_sum .*= Float32(365.25)   # in cm/year
    path_out = joinpath(dirout, "sbas_rate.un")
    write_unw_igram!(raw, path_out, amp_sum, psidt_sum)
end

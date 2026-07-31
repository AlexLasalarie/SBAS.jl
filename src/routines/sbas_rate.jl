"""
    sbas_rate(
        ts_list::String,
        wid::Int,
        len::Int
    )

Compute the linear surface deformation rate (in cm/year) from the SBAS time 
series. 

# Positional Arguments
- `ts_list`: `String` path to file containing the list of deformation files
- `wid`: `Int` width of unwrapped interferogram files (in pixels)
- `len`: `Int` length of unwrapped interferogram files (in pixels)

# Returns
`sbas_rate.un`: flat binary of `.un` format containing the `RMSE` in the first 
band and the `linear rate` (in cm/year) in the second band.
"""
function sbas_rate(
    ts_list::String,
    wid::Int,
    len::Int
)
    # Read the list of files
    files = parse_file(ts_list, "s")[1]
    npts = length(files)

    # Extract dates and sort list
    dates = Vector{Date}(undef, npts)
    dt = Vector{Float32}(undef, npts)
    for (k, file) in enumerate(files)
        date1, date2 = extract_dates(file)
        dates[k] = date2
        dt[k] = Float32((date2 - date1).value)
    end

    # Buffers
    raw = Matrix{Float32}(undef, 2 * wid, len)
    amp = Matrix{Float32}(undef, len, wid)
    def = Matrix{Float32}(undef, len, wid)
    cube = Matrix{Float32}(undef, npts, len * wid)
    def_rate = Matrix{Float32}(undef, len, wid)
    rmse = Matrix{Float32}(undef, len, wid)

    # Fill cube
    for (k, file) in enumerate(files)
        read_unw_igram!(raw, amp, def, file)
        for idx in eachindex(def)
            cube[k, idx] = def[idx]
        end
    end

    # LS inversion - y=mx+b
    sumx = sum(dt)
    sumx2 = dot(dt, dt)
    den = npts * sumx2 - sumx * sumx
    wvec = Vector{Float32}(undef, npts)
    for (k, xk) in enumerate(dt)
        wvec[k] = npts * xk - sumx
    end

    # Compute deformation rate
    @threads for idx in 1:len*wid

        # Slope
        yvec = view(cube, :, idx)
        num = Float32(0.0)
        sumy = Float32(0.0)
        for (yk, wk) in zip(yvec, wvec)
            num += wk * yk
            sumy += yk
        end
        m = num / den                       # slope [cm/day]

        # Intercept
        b = (sumy - m * sumx) / npts        # intercept [cm]

        # RMSE
        sum_e2 = Float32(0.0)
        for (xk, yk) in zip(dt, yvec)
            fit = m * xk + b
            err = fit - yk
            sum_e2 += err * err
        end

        # Populate
        def_rate[idx] = m * 365.25          # deformation rate [cm/year]
        rmse[idx] = sqrt(sum_e2 / npts)     # RMSE [cm]
    end

    # Output the deformation rate
    dirout = dirname(ts_list)
    pathout = joinpath(dirout, "sbas_rate.un")
    write_unw_igram!(raw, pathout, rmse, def_rate)
    return rmse, def_rate
end

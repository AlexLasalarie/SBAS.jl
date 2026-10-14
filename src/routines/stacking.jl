"""
    stacking(
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
function stacking(
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

"""
Stacking with outlier rejection

Common scene (k) sum over all pairs (n):
    v_k = sum(|d_kn|) / sum(|dt_kn|)
"""
function stacking_mad(
    wid::Int,
    len::Int,
    ref_col::Int,
    ref_row::Int,
    intlist::String;
    wavelength::Real=5.55
)
    # Read list of files
    ints, data = metadata(wid, len, intlist)

    # Buffers
    raw = Matrix{Float32}(undef, 2 * wid, len)
    amp = Matrix{Float32}(undef, len, wid)
    psi = Matrix{Float32}(undef, len, wid)
    vel_buf = Array{Float32}(undef, len, wid, data.nacq)
    abs_dev = Vector{Float32}(undef, data.nacq)
    isgood = BitArray(undef, len, wid, data.nacq)
    amp_sum = zeros(Float32, len, wid)
    psi_sum = zeros(Float32, len, wid)
    psidt_sum = zeros(Float32, len, wid)
    dtsum_day = zeros(Float32, len, wid)
    dt2sum_day = zeros(Float32, len, wid)

    # Common scene velocities
    for (k, scene) in enumerate(data.scenes)
        # Reset
        fill!(psi_sum, zero(Float32))
        dt_sum = zero(Int)
        for int in ints
            if (int.date1 == scene) || (int.date2 == scene)
                read_unw_igram!(raw, amp, psi, int.name)
                calibrate!(psi, ref_row, ref_col)
                psi_sum .+= abs.(psi)
                dt_sum += abs(int.dt)
            end
        end
        vel = view(vel_buf, :, :, k)
        vel .= psi_sum ./ dt_sum
    end

    # MAD
    for col in 1:wid
        for row in 1:len
            vel_vec = view(vel_buf, row, col, :)
            med = median(vel_vec)
            for (k, vel) in enumerate(vel_vec)
                abs_dev[k] = abs(vel - med)
            end
            mad = median(abs_dev)
            # max_vel = med + 4 * 1.483 * mad         # proper normalization
            max_vel = med + 4 * mad                 # in Ke's code
            for (k, vel) in enumerate(vel_vec)
                isgood[row, col, k] = vel <= max_vel
            end
        end
    end

    # Stacking with a switch
    fill!(psi_sum, zero(Float32))
    for int in ints

        # Read interferogram
        read_unw_igram!(raw, amp, psi, int.name)
        idx1 = int.idx1
        idx2 = int.idx2
        dt = int.dt
        if int.flip
            psi .*= Float32(-1.0)
            dt = abs(dt)
        end

        # Add to running sums if valid (check at each pixel)
        amp_sum .+= amp
        for col in 1:wid
            for row in 1:len
                if isgood[row, col, idx1] && isgood[row, col, idx2]
                    psi_val = psi[row, col]
                    psi_sum[row, col] += psi_val
                    psidt_sum[row, col] += psi_val * dt
                    dtsum_day[row, col] += dt
                    dt2sum_day[row, col] += dt * dt
                end
            end
        end
    end
    amp_sum ./= data.nint

    # Stacking deformation rate [cm/year]
    coef = Float32(wavelength / (4 * pi))
    psi_sum .*= coef                        # [cm]
    psi_sum ./= dtsum_day                   # [cm/day]
    psi_sum .*= Float32(365.25)             # [cm/year]
    dirout = dirname(ints[1].name)
    path_out = joinpath(dirout, "stack_rate_mad.un")
    write_unw_igram!(raw, path_out, amp_sum, psi_sum)

    # Linear SBAS deformation rate [cm/year]
    psidt_sum .*= coef                      # [cm * day]
    psidt_sum ./= dt2sum_day                # [cm/day]
    psidt_sum .*= Float32(365.25)           # [cm/year]
    path_out = joinpath(dirout, "sbas_rate_mad.un")
    write_unw_igram!(raw, path_out, amp_sum, psidt_sum)
end

"""
    sbas_list(
        max_err::Real,
        max_atm::Real,
        path_err::String,
        path_atm::String;
        name_out::String="sbas_pairs"
    )

Create the list of interferograms to use for SBAS inversion by imposing a 
baseline on unwrapping error and atmospheric noise.

# Positional Arguments
- `max_err`: `Real` maximum unwrapping error (in rad)
- `max_atm`: `Real` maximum atmospheric delay (in rad)
- `path_err`: `String` path to file containing the unwrapping errors
- `path_atm`: `String` path to file containing the atmospheric delays

# Optional Arguments
- `name_out`: `String` path to output file (default is `"sbas_pairs"`)

# Returns
List of unwrapped interferograms to use for SBAS inversion and path to the list.
"""
function sbas_list(
    max_err::Real,
    max_atm::Real,
    path_err::String,
    path_atm::String;
    name_out::String="sbas_pairs"
)
    # Read lists
    path_ints, unw_errs = parse_file(path_err, "sf")
    path_acqs, atm_errs = parse_file(path_atm, "sf")
    nint = length(path_ints)
    nacq = length(path_acqs)

    # Statistics
    qerr = round.(Int, quantile(unw_errs, (0.0, 0.25, 0.5, 0.75, 1.0)))
    qatm = round.(Int, quantile(atm_errs, (0.0, 0.25, 0.5, 0.75, 1.0)))

    # Determine if scene is valid
    isvalid_acq = Dict{Date,Bool}()
    nval_scenes = 0
    for (path, atm) in zip(path_acqs, atm_errs)
        date = extract_dates(path)[1]
        isvalid_acq[date] = (atm < max_atm)
        if atm < max_atm
            nval_scenes += 1
        end
    end

    # Determine if interferogram is valid
    used_scenes = Set{Date}()
    adj_list = Dict{Date,Set{Date}}()
    isvalid_int = Vector{Bool}(undef, nint)
    nval_ints = 0
    nval_pairs = 0
    for (k, path_int) in enumerate(path_ints)
        dates = extract_dates(path_int)
        isvalid_date1 = get(isvalid_acq, dates[1], false)
        isvalid_date2 = get(isvalid_acq, dates[2], false)
        isvalid = (unw_errs[k] < max_err) && isvalid_date1 && isvalid_date2
        isvalid_int[k] = isvalid
        if isvalid
            nval_pairs += 1
            push!(used_scenes, dates[1])
            push!(used_scenes, dates[2])
            push!(get!(adj_list, dates[1], Set{Date}()), dates[2])
            push!(get!(adj_list, dates[2], Set{Date}()), dates[1])
        end
        if unw_errs[k] < max_err
            nval_ints += 1
        end
    end
    isconnected = network_connection(used_scenes, adj_list)

    # Output the statistics
    dirout = dirname(path_err)
    path_out = joinpath(dirout, "partition_summary")
    open(path_out, "w") do io
        println(io, "# Unwrapping error")
        println(io, "min, Q1, Q2, Q3, max: $(qerr[1]), $(qerr[2]), $(qerr[3]), $(qerr[4]), $(qerr[5])")
        println(io, "Unwrapping error baseline: $(max_err)")
        println(io, "Valid interferograms: $(nval_ints) of $(nint)\n")
        println(io, "# Atmospheric delays")
        println(io, "min, Q1, Q2, Q3, max: $(qatm[1]), $(qatm[2]), $(qatm[3]), $(qatm[4]), $(qatm[5])")
        println(io, "Maximum atm. standard dev.: $(max_atm)")
        println(io, "Valid scenes: $(nval_scenes) of $(nacq)\n")
        println(io, "# Used for deformation")
        println(io, "Number of scenes: $(length(used_scenes))")
        println(io, "Number of pairs: $(nval_pairs) of $(nint)")
        println(io, "Is network connected? $(isconnected)")
    end

    # Output to list
    path_out = joinpath(dirout, name_out)
    open(path_out, "w") do io
        for (k, path) in enumerate(path_ints)
            if isvalid_int[k]
                println(io, path)
            end
        end
    end
    return path_out
end

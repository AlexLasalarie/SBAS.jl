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
- `path_out`: `String` path to output file (default is `"sbas_pairs"`)

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

    # Determine if scene is valid
    isvalid_acq = Dict{Date,Bool}()
    for (path, atm) in zip(path_acqs, atm_errs)
        date = extract_dates(path)[1]
        isvalid_acq[date] = (atm < max_atm)
    end

    # Determine if interferogram is valid
    isvalid_int = Vector{Bool}(undef, nint)
    for (k, path_int) in enumerate(path_ints)
        dates = extract_dates(path_int)
        isvalid_date1 = get(isvalid_acq, dates[1], false)
        isvalid_date2 = get(isvalid_acq, dates[2], false)
        isvalid_int[k] = (unw_errs[k] < max_err) && isvalid_date1 && isvalid_date2
    end

    # Apply baselines
    dirout = dirname(path_err)
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

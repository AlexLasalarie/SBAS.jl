function extract_dates(path::String)
    pattern = r"(\d{8})(?:_(\d{8}))?\.[a-zA-Z0-9]+$"
    m = match(pattern, path)
    if isnothing(m)
        return Date[]
    end
    date_format = dateformat"yyyymmdd"
    dates = [Date(d, date_format) for d in m.captures if !isnothing(d)]
    return dates
end

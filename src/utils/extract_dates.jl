"""
    extract_dates(path::String)

Extract `date1` and `date2` from the string `/path/to/date1_date2.ext`.

Only the basename is searched, so dates appearing in parent directories are
ignored. Any prefix and any number of extensions are tolerated, as in
`filt_date1_date2.int.unw`. When more than two dates remain, the two nearest the 
extension are taken.

The dates are returned in the order they appear in the name, *not* in
chronological order: a pair named `date2_date1` is reported reversed, and it is
the caller's job to handle the sign (see `metadata`).

# Positional Arguments
- `path`: `String` path to the file

# Returns
`Vector{Date}` holding the pair of dates, a single date for names of the form
`date.ext`, or an empty vector when the name carries no date at all.
"""
function extract_dates(path::String)
    date_format = dateformat"yyyymmdd"
    dates = Date[]
    for m in eachmatch(r"(?<![0-9])([0-9]{8})(?![0-9])", basename(path))
        date = tryparse(Date, m[1], date_format)
        isnothing(date) || push!(dates, date)
    end
    return length(dates) > 2 ? dates[end-1:end] : dates
end

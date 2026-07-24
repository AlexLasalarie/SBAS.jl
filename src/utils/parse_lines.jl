function parse_file(
    filename::String,
    format_spec::String;
    delim=r"\s+")
    type_map = Dict(
        's' => (String, s -> String(s)),
        'f' => (Float64, s -> parse(Float64, s)),
        'i' => (Int, s -> parse(Int, s))
    )
    cols = Tuple([type_map[c][1][] for c in format_spec])
    parsers = [type_map[c][2] for c in format_spec]
    open(filename, "r") do io
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(strip(line), delim)
            for (i, parser) in enumerate(parsers)
                push!(cols[i], parser(parts[i]))
            end
        end
    end
    return cols
end

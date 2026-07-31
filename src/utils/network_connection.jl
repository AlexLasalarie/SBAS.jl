function network_connection(
    used_scenes::Set{Date},
    adj_list::Dict{Date,Set{Date}}
)::Bool
    if isempty(used_scenes)
        return true
    end
    visited = Set{Date}()
    start_node = first(used_scenes)
    queue = [start_node]
    push!(visited, start_node)
    while !isempty(queue)
        curr = popfirst!(queue)
        for neighbor in get(adj_list, curr, Set{Date}())
            if !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    return length(visited) == length(used_scenes)
end

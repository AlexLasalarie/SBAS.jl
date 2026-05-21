struct IgramInfo
    name::String
    wid::Int
    len::Int
    date1::Date
    date2::Date
    idx1::Int
    idx2::Int
end

struct DataInfo
    wid::Int
    len::Int
    npix::Int
    nint::Int
    nacq::Int
    dir::String
    scenes::Vector{Date}
    out_paths::Vector{String}
end

struct SystemInfo
    total_ram::Float64
    avail_ram::Float64
    threads::Int
end

struct ChunkInfo
    path::String
    p1::Int
    p2::Int
    numel::Int
end

struct TileInfo
    path::String
    offset::Vector{Int}
end

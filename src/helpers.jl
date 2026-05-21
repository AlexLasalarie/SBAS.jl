# Compile information on each interferogram and on the overall data stack
function metadata(
    wid::Int,
    len::Int,
    intlist::String
)
    # Read the file list
    file_paths = Vector{String}()
    open(intlist, "r") do io
        for line in eachline(io)
            if !isempty(strip(line))
                # Check the file exists
                if isfile(line)
                    push!(file_paths, line)
                else
                    error("Could not find: $(line)")
                end
            end
        end
    end
    nint = length(file_paths)

    # Form list of acquisition times 
    d1 = Vector{String}(undef, nint)
    d2 = Vector{String}(undef, nint)
    pattern = r"(\d{8})_(\d{8})"
    for (i, file_path) in enumerate(file_paths)
        m = match(pattern, file_path)
        d1[i] = String(m[1])
        d2[i] = String(m[2])
    end

    # Unique scenes
    unique_scenes = sort(unique([d1; d2]))
    nacq = length(unique_scenes)
    scenes = Date.(unique_scenes, "yyyymmdd")

    # Form list of indices
    idx1 = Vector{Int}(undef, nint)
    idx2 = Vector{Int}(undef, nint)
    for i in 1:nint
        idx1[i] = findfirst(x -> x == d1[i], unique_scenes)
        idx2[i] = findfirst(x -> x == d2[i], unique_scenes)
    end

    # Bundle dataset information
    dir_in = dirname(abspath(file_paths[1]))
    dir_out = joinpath(dir_in, "time_series")
    ref_scene = unique_scenes[1]
    out_paths = Vector{String}(undef, nacq - 1)
    if isdir(dir_out)
        rm(dir_out, recursive=true)
    end
    mkdir(dir_out)
    for (i, scene) in enumerate(unique_scenes[2:end])
        out_paths[i] = joinpath(dir_out, "sbas_$(ref_scene)_$(scene).un")
    end
    data = DataInfo(wid, len, wid * len, nint, nacq, dir_in, scenes, out_paths)

    # Create object containing the igram information
    ints = Vector{IgramInfo}(undef, nint)
    for i in 1:nint
        ints[i] = IgramInfo(
            file_paths[i],
            wid,
            len,
            Date.(d1[i], "yyyymmdd"),
            Date.(d2[i], "yyyymmdd"),
            idx1[i],
            idx2[i]
        )
    end
    return ints, data
end

# Generate the B^-1 matrix
function make_design_matrix(ints::Vector{IgramInfo}, data::DataInfo)
    B = zeros(Float64, data.nint, data.nacq - 1)
    tkp1mtk = data.scenes[2:end] .- data.scenes[1:end-1]
    days_f64 = Float64.(Dates.value.(tkp1mtk))
    for i in 1:data.nint
        col1 = ints[i].idx1
        col2 = ints[i].idx2
        for j in col1:col2-1
            B[i, j] = days_f64[j]
        end
    end
    Binv = pinv(B)
    Binv_f32 = Float32.(Binv)
    return Binv_f32, Float32.(days_f64)
end

# Collect information on system 
function system_info(fraction_allowed)
    total_ram = Float64(Sys.total_memory())
    n_threads = Threads.nthreads()
    avail_ram = total_ram * fraction_allowed
    sys = SystemInfo(total_ram, avail_ram, n_threads)
    return sys
end

# Handles the batching logic
function batch_logic(
    data::DataInfo,
    sys::SystemInfo
)

    # Chunking logic -- each chunk contains the full insar measurement history at
    # a subset of the pixels. Each chunk is numel x nint where numel is the
    # number of elements (pixels) in the chunk
    pix_per_chunk = max(1, floor(Int, sys.avail_ram / (4 * (data.nint + data.nacq - 1))))
    nchunk = ceil(Int, data.npix / pix_per_chunk)
    chunk_dir = joinpath(data.dir, "chunks")
    if isdir(chunk_dir)
        rm(chunk_dir, recursive=true)
    end
    mkdir(chunk_dir)
    chunks = Vector{ChunkInfo}(undef, nchunk)
    for i in 1:nchunk
        p1 = pix_per_chunk * (i - 1) + 1                # start pixel index
        p2 = min(i * pix_per_chunk, data.npix)          # end pixel index
        chunk_path = joinpath(chunk_dir, "chunk_$i")
        touch(chunk_path)
        numel = p2 - p1 + 1
        chunks[i] = ChunkInfo(chunk_path, p1, p2, numel)
    end

    # Tiling logic -- the tile contains the time series at all pixels. The tile
    # is npix x (nacq-1), where nacq is the number of acquisitions.
    tile_path = abspath(joinpath(chunk_dir, "time_series"))
    offsets = Vector{Int}(undef, data.nacq - 1)
    for i in 1:data.nacq-1
        offsets[i] = (i - 1) * data.npix
    end
    tile = TileInfo(tile_path, offsets)
    touch(tile_path)

    return chunks, tile
end

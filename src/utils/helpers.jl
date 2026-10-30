# Compile information on each interferogram and on the overall data stack
function metadata(
    wid::Int,
    len::Int,
    intlist::String
)
    # Read the file list
    @assert isfile(intlist) "$intlist is not a valid file"
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
    d1 = Vector{Date}(undef, nint)
    d2 = Vector{Date}(undef, nint)
    for (i, file_path) in enumerate(file_paths)
        dates = extract_dates(file_path)
        if length(dates) != 2
            error("Could not extract a pair of dates from: $(file_path)")
        end
        d1[i] = dates[1]
        d2[i] = dates[2]
    end

    # Unique scenes (sorted chronologically)
    scenes = sort(unique([d1; d2]))
    nacq = length(scenes)

    # Form list of indices
    idx1 = Vector{Int}(undef, nint)
    idx2 = Vector{Int}(undef, nint)
    for i in 1:nint
        idx1[i] = findfirst(==(d1[i]), scenes)
        idx2[i] = findfirst(==(d2[i]), scenes)
    end

    # Bundle dataset information
    dir_in = dirname(abspath(file_paths[1]))
    dir_out = joinpath(dir_in, "time_series")
    ref_scene = Dates.format(scenes[1], "yyyymmdd")
    out_paths = Vector{String}(undef, nacq - 1)
    if isdir(dir_out)
        rm(dir_out, recursive=true)
    end
    mkdir(dir_out)
    for (i, scene) in enumerate(scenes[2:end])
        scene_str = Dates.format(scene, "yyyymmdd")
        out_paths[i] = joinpath(dir_out, "sbas_$(ref_scene)_$(scene_str).un")
    end
    data = DataInfo(wid, len, wid * len, nint, nacq, dir_in, scenes, out_paths)

    # Create object containing the igram information
    ints = Vector{IgramInfo}(undef, nint)
    for i in 1:nint
        date1 = d1[i]
        date2 = d2[i]
        col1 = idx1[i]
        col2 = idx2[i]

        # Normalize to chronological order. Pairs named date2_date1 carry the
        # opposite sign, so flag them and let the readers negate the phase.
        flip = date2 < date1
        if flip
            date1, date2 = date2, date1
            col1, col2 = col2, col1
        end
        dt = Int((date2 - date1).value)
        ints[i] = IgramInfo(
            file_paths[i],
            wid,
            len,
            date1,
            date2,
            col1,
            col2,
            dt,
            flip
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

    # Flag disconnected network
    nparam = data.nacq - 1
    rankB = rank(B)
    if rankB < nparam
        @warn "Design matrix is rank deficient: the interferogram network does " *
              "not constrain every epoch. The time series will be damped toward " *
              "zero across the unconstrained intervals." rank = rankB parameters = nparam
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

"""
    common_scene_stack(
        wid::Int,
        len::Int,
        intlist::String;
        dtmax::Int=60,
        name_out::String="atmospheric_delays"
    )

Estimate atmospheric contributions using common scene stacking.

# Positional Arguments
- `wid`: `Int` width of the unwrapped interfergorams (in pixels)
- `len`: `Int` length of the unwrapped interferograms (in pixels)
- `intlist`: `String` path to the file containing the list of files

# Optional arguments
- `dtmax`: `Int` maximum temporal baseline for stacking (default is 60 days)
- `name_out`: `String` name of output file (default is "atmospheric_delays"

# Returns
1. List of the standard deviation of the common scene stacks
2. Write common scene stacks to disk as unwrapped interferogram
"""
function common_scene_stack(
    wid::Int,
    len::Int,
    intlist::String,
    dtmax::Int=60;
    name_out::String="atmospheric_delays"
)
    # Pre-allocate buffers
    raw = Matrix{Float32}(undef, 2 * wid, len)
    amp = Matrix{Float32}(undef, len, wid)
    psi = Matrix{Float32}(undef, len, wid)
    css = Matrix{Float32}(undef, len, wid)
    amp_sum = Matrix{Float32}(undef, len, wid)

    # Read list of files
    ints, data = metadata(wid, len, intlist)

    # Create the output directory
    outputdir = joinpath(data.dir, "common_scene_stacks")
    if isdir(outputdir)
        rm(outputdir; recursive=true)
    end
    mkdir(outputdir)

    # Compute common scene stack
    println("Computing stacks...")
    css_list = Vector{Float32}(undef, data.nacq)
    paths = Vector{String}(undef, data.nacq)
    for (k, scene) in enumerate(data.scenes)

        # Reset
        fill!(css, Float32(0.0))
        fill!(amp_sum, Float32(0.0))

        # Forward sum
        fcnt = 0
        for int in ints
            if (int.idx2 != k) || (int.dt > dtmax)
                continue
            end
            fcnt += 1
            read_unw_igram!(raw, amp, psi, int.name)
            css .-= psi
            amp_sum .+= amp
        end

        # Backward sum
        bcnt = 0
        for int in ints
            if (int.idx1 != k) || (int.dt > dtmax)
                continue
            end
            bcnt += 1
            read_unw_igram!(raw, amp, psi, int.name)
            css .+= psi
            amp_sum .+= amp
        end

        # Result
        date_str = Dates.format(scene, "yyyymmdd")
        paths[k] = joinpath(outputdir, "$(date_str).un")
        css ./= (fcnt + bcnt)
        amp_sum ./= (fcnt + bcnt)
        write_unw_igram!(raw, paths[k], amp_sum, css)
        css_list[k] = std(css)
    end

    # Write to file
    println("Writing to summary file...")
    path_out = joinpath(data.dir, name_out)
    open(path_out, "w") do io
        for (k, val) in enumerate(css_list)
            println(io, "$(paths[k]) $(val)")
        end
    end
end

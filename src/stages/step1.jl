function step1(
    data::DataInfo,
    ints::Vector{IgramInfo},
    chunks::Vector{ChunkInfo},
    ref_row::Int,
    ref_col::Int
)
    # Greet
    println("Forming $(length(chunks)) temporary chunks...")

    # Reusable buffers
    unw_buffer = Matrix{Float32}(undef, data.len, data.wid)         # Phase buffer 
    raw_buffer = Matrix{Float32}(undef, 2 * data.wid, data.len)     # row-major read buffer 
    amp_buffer = Matrix{Float32}(undef, data.len, data.wid)         # amplitude buffer
    amp_cumsum = zeros(Float32, data.len, data.wid)                 # for average amplitude

    # Open all chunk files
    ios = [open(chunk.path, "w") for chunk in chunks]

    # Sequential read igram, write to chunk
    try
        for i in 1:data.nint

            # Read a single igram
            read_unw_igram!(raw_buffer, amp_buffer, unw_buffer, ints[i].name)
            calibrate!(unw_buffer, ref_row, ref_col)

            # Write to chunks
            for (k, chunk) in enumerate(chunks)
                chunk_slice = view(unw_buffer, chunk.p1:chunk.p2)
                write(ios[k], chunk_slice)
            end

            # Update amplitude cumsum 
            amp_cumsum .+= amp_buffer

        end
    finally
        foreach(close, ios)
    end
    println("Formed the temporary chunks.")
    return amp_cumsum ./ data.nint
end



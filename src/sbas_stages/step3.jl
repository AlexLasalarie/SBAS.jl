function step3(avg_amp::Matrix{Float32}, tile::TileInfo, data::DataInfo, ref_row::Int, ref_col::Int)

    # Greet 
    println("Forming deformation maps")

    # Buffers 
    raw_buffer = Matrix{Float32}(undef, 2 * data.wid, data.len)
    phase = Matrix{Float32}(undef, data.len, data.wid)

    # The read is done little by little
    open(tile.path, "r") do io
        for out_path in data.out_paths
            read!(io, phase)
            calibrate!(phase, ref_row, ref_col)
            write_unw_igram!(raw_buffer, out_path, avg_amp, phase)
        end
    end

    # Clean up by removing the temporary folder
    chunk_dir = dirname(tile.path)
    if isdir(chunk_dir)
        rm(chunk_dir, recursive=true)
    end
end


"""
    read_unw_igram!(
        raw_buffer::AbstractMatrix{Float32}, 
        amp_buffer::AbstractMatrix{Float32}, 
        phase_buffer::AbstractMatrix{Float32}, 
        filename::String)

Read an unwrapped interferogram amplitude and phase directly into a pre-allocated buffer.

# Arguments
- `raw_buffer`: pre-allocated `AbstractMatrix{Float32}` of size `(2*width, length)`.
- `amp_buffer`: pre-allocated `AbstractMatrix{Float32}` of size `(length, width)`.
- `phase_buffer`: pre-allocated `AbstractMatrix{Float32}` of size `(length, width)`.
- `filename`: path to the binary interferogram file.

# Returns
- `raw_buffer`: `Matrix{Float32}` raw data (overwrites the buffer).
- `amp_buffer`: `Matrix{Float32}` amplitude of interferogram (overwrites the buffer).
- `phase_buffer`: `Matrix{Float32}` phase of interferogram (overwrites the buffer).

# Note
This function is the high-performance, in-place version of `read_unw_igram`.
Use this version for batch processing where the buffer can be reused.
The input buffers will be overwritten!
"""
function read_unw_igram!(
    raw_buffer::AbstractMatrix{Float32},
    amp_buffer::AbstractMatrix{Float32},
    phase_buffer::AbstractMatrix{Float32},
    file_name::String
)

    # Bound check 1
    len, wid = size(amp_buffer)
    len1, wid1 = size(phase_buffer)
    if (len != len1) || (wid != wid1)
        error("The amplitude and phase matrices are not the same size!")
    end

    # Bound check 2 (raw_buffer is size (2*width, length))
    len2, wid2 = size(raw_buffer)
    if (len != wid2) || (2 * wid != len2)
        error("The raw buffer is not the correct size (must be (2*width, length))!")
    end

    # Fill raw_buffer
    open(file_name, "r") do io
        read!(io, raw_buffer)
    end

    # Fill amp and phase (no allocations)
    _, wid = size(amp_buffer)
    amp_view = view(raw_buffer, 1:wid, :)
    phase_view = view(raw_buffer, (wid+1):(2*wid), :)
    amp_buffer .= amp_view'
    phase_buffer .= phase_view'

    return nothing
end



"""
    amp, phase = read_unw_igram(filename::String, wid::Int, len::Int)

Read an unwrapped interferogram amplitude and phase.

# Arguments
- `filename`: Path to wrapped interferogram file.
- `igram_width`: `Int` width of the input interferogram.
- `igram_length`: `Int` length of the input interferogram.

# Returns
- `amp`: `Matrix{Float32}` amplitude of unwrapped interferogram. 
- `phase`: `Matrix{Float32}` phase of unwrapped interferogram. 
"""
function read_unw_igram(file_name, igram_width, igram_length)
    raw_buffer = Matrix{Float32}(undef, 2 * igram_width, igram_length)
    amp_buffer = Matrix{Float32}(undef, igram_length, igram_width)
    phase_buffer = Matrix{Float32}(undef, igram_length, igram_width)
    read_unw_igram!(raw_buffer, amp_buffer, phase_buffer, file_name)
    return amp_buffer, phase_buffer
end

"""
    write_unw_igram!(
        raw_buffer::AbstractMatrix{Float32}, 
        name_out::String, 
        amp::AbstractMatrix{Float32}, 
        phase::AbstractMatrix{Float32}
    )

Write the amplitude and phase of an unwrapped interferogram to a binary file.
The data is stacked as [amplitude'; phase'].

# Arguments
- `raw_buffer`: pre-allocated `AbstractMatrix{Float32}` of size `(2*width, length)`.
- `name_out`: path to the binary interferogram file.
- `amp`: `AbstractMatrix{Float32}` of size `(length, width)`.
- `phase`: `AbstractMatrix{Float32}` of size `(length, width)`.

# Returns
- `AbstractMatrix{Float32}`: modified `raw_buffer`.

# Note
This function is the high-performance, in-place version of `write_unw_igram`.
Use this version for batch processing where the buffer can be reused.
The raw_buffer must be of size (2*width, length)!
"""
function write_unw_igram!(
    raw_buffer::AbstractMatrix{Float32},
    name_out::String,
    amp::AbstractMatrix{Float32},
    phase::AbstractMatrix{Float32}
)

    # Bound check 1
    len, wid = size(amp)
    len1, wid1 = size(phase)
    if (len != len1) || (wid != wid1)
        error("The amplitude and phase matrices are not the same size!")
    end

    # Bound check 2 (raw_buffer is size (2*width, length))
    len2, wid2 = size(raw_buffer)
    if (len != wid2) || (2 * wid != len2)
        error("The raw buffer is not the correct size (must be (2*width, length))!")
    end

    # Fills buffer (no allocs)
    amp_view = view(raw_buffer, 1:wid, :)
    phase_view = view(raw_buffer, (wid+1):(2*wid), :)
    amp_view .= amp'
    phase_view .= phase'

    # Writes to disk
    open(name_out, "w") do io
        write(io, raw_buffer)
    end

    # Exit
    return nothing
end

"""
    write_unw_igram(name_out::String, amp::AbstractMatrix{Float32}, phase::AbstractMatrix{Float32})

Write the amplitude and phase of an unwrapped interferogram to a binary file.
The data is stacked as [amplitude'; phase'].

# Arguments
- `name_out`: Path to unwrapped interferogram file.
- `amp`: `AbstractMatrix{Float32}` of size `(length, width)`.
- `phase`: `AbstractMatrix{Float32}` of size `(length, width)`.
"""
function write_unw_igram(name_out::String, amp::AbstractMatrix{Float32}, phase::AbstractMatrix{Float32})

    # Creates buffer
    len, wid = size(amp)
    raw_buffer = Matrix{Float32}(undef, 2 * wid, len)

    # Writes to file
    write_unw_igram!(raw_buffer, name_out, amp, phase)
end

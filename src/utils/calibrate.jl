"""
    calibrate!(
        phase::AbstractMatrix{Float32}, 
        ref_row::Int, 
        ref_col::Int; 
        box_size::Int=3
    )

Calibrate an unwrapped inteferogram to a reference location. Note that the `phase` 
will be overwritten in RAM.

# Arguments
- `phase`: `AbstractMatrix{Float32}` phase of unwrapped interferogram (will be overwritten)
- `ref_row`: `Int` row location of reference pixel
- `ref_col`: `Int` column location of reference pixel 
- `box_size`: `Int` size of box around reference pixel (default is `3` for a 3x3 box)

# Returns
- `Matrix{Float32}`: overwrites the initial matrix with calibrated phase matrix

# Note
The input matrix gets overwritten in RAM!
"""
function calibrate!(
    phase::AbstractMatrix{Float32},
    ref_row::Int,
    ref_col::Int;
    box_size::Int=3
)

    # Dimensions
    len, wid = size(phase)

    # Reference cluster
    del = floor(Int, box_size / 2)
    row1 = max(ref_row - del, 1)
    row2 = min(ref_row + del, len)
    col1 = max(ref_col - del, 1)
    col2 = min(ref_col + del, wid)
    npix = (row2 - row1 + 1) * (col2 - col1 + 1)

    # Mean of cluster
    avg = 0.f0
    for j in col1:col2
        for i in row1:row2
            avg += phase[i, j]
        end
    end

    # Calibrates
    avg /= npix
    phase .-= avg
end



"""
    calibrate!(signal::AbstractMatrix{ComplexF32}, ref_row::Int, ref_col::Int; box_size::Int=3)

Calibrate a wrapped inteferogram to a reference location. Note that the `signal` matrix 
will be overwritten in RAM.

# Arguments
- `signal`: `AbstractMatrix{ComplexF32}` wrapped interferogram (will be overwritten)
- `ref_row`: `Int` row location of reference pixel
- `ref_col`: `Int` column location of reference pixel 
- `box_size`: `Int` size of box around reference pixel (default is `3` for a 3x3 box)

# Returns
- `Matrix{ComplexF32}`: overwrites the initial wrapped interferogram with calibrated version

# Note
The input matrix gets overwritten in RAM!
"""
function calibrate!(signal::AbstractMatrix{ComplexF32}, ref_row::Int, ref_col::Int; box_size::Int=3)

    # Dimensions
    len, wid = size(signal)

    # Reference cluster
    del = floor(Int, box_size / 2)
    row1 = max(ref_row - del, 1)
    row2 = min(ref_row + del, len)
    col1 = max(ref_col - del, 1)
    col2 = min(ref_col + del, wid)

    # Mean of cluster
    phasor = ComplexF32(0.0)
    for j in col1:col2
        for i in row1:row2
            s_ij = signal[i, j]
            mag_sij = abs(s_ij)
            if mag_sij > 0
                phasor += s_ij / mag_sij
            end
        end
    end

    # Calibrates
    if abs(phasor) > 0
        phasor /= abs(phasor)
    else
        phasor = ComplexF32(1.0)
    end
    signal .*= conj(phasor)
end

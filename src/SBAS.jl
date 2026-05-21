module SBAS

# External dependencies
using LinearAlgebra
using Base.Threads
using Dates

# Internal dependencies

# Types
include("types.jl")

# File IO
include("fileio/read_unw_igram.jl")
include("fileio/write_unw_igram.jl")

# Helpers
include("helpers.jl")
include("calibrate.jl")

# Stages
include("stages/step1.jl")
include("stages/step2.jl")
include("stages/step3.jl")

# Exports
export sbas

# Body
"""
    sbas()
"""
function sbas(
    wid::Int,
    len::Int,
    ref_row::Int,
    ref_col::Int,
    intlist::String;
    fraction_ram::Float64=0.1,
    wavelength::Float64=5.55
)

    # Compile information about the dataset and individual interferograms
    ints, data = metadata(wid, len, intlist)

    # Form the design matrix
    Binv, tkp1mtk = make_design_matrix(ints, data)

    # System information
    sys = system_info(fraction_ram)

    # Batching logic
    chunks, tile = batch_logic(data, sys)

    # Read the files sequentially and populate chunks
    avg_amp = step1(data, ints, chunks)

    # Solve for the time series
    step2(data, chunks, tile, sys, Binv, wavelength, tkp1mtk)

    # Form the maps
    step3(avg_amp, tile, data, ref_row, ref_col)
end

end

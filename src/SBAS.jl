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
    sbas(
        wid::Int,
        len::Int,
        ref_row::Int,
        ref_col::Int,
        intlist::String;
        fraction_ram::Float64=0.1,
        wavelength::Float64=5.55
    )

Generate an InSAR deformation time series from a stack of unwrapped interferograms
using the SBAS algorithm.

# Positional Arguments
- `wid`: `Int` width of the interferograms (in pixels)
- `len`: `Int` length of the interferograms (in pixels)
- `ref_row`: `Int` row position of reference pixel (in pixels)
- `ref_col`: `Int` column position of referece pixel (in pixels)
- `intlist`: `String` path to the file containing the list of unwrapped interferograms

# Optional Arguments 
- `fraction_ram`: `Float64` fraction of random access memory allocated to the run 
    (default is 0.1 or 10% of full available RAM)
- `wavelength`: `Float64` wavelength of sensor, in cm (default is C-Band: 5.5 cm)

# Returns
Deformation time series. If the network of interferograms has n scenes, n-1 maps of
the deformation (in cm) with respect to the first acquistion will be generated. 

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
    avg_amp = step1(data, ints, chunks, ref_row, ref_col)

    # Solve for the time series
    step2(data, chunks, tile, sys, Binv, wavelength, tkp1mtk)

    # Form the maps
    step3(avg_amp, tile, data, ref_row, ref_col)

    println("Done!")
end

end

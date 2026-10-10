module SBAS

# External dependencies
using LinearAlgebra
using Base.Threads
using Dates
using Statistics        # for common scene stacking

# Internal dependencies

# Types
include("types.jl")

# File IO
include("fileio/read_unw_igram.jl")
include("fileio/write_unw_igram.jl")

# Utils
include("utils/parse_lines.jl")
include("utils/extract_dates.jl")
include("utils/calibrate.jl")
include("utils/helpers.jl")
include("utils/network_connection.jl")

# Routines
include("routines/mask.jl")
include("routines/unwrapping_error.jl")
include("routines/common_scene_stack.jl")
include("routines/sbas_list.jl")
include("routines/stacking.jl")

# SBAS stages
include("sbas_stages/step1.jl")
include("sbas_stages/step2.jl")
include("sbas_stages/step3.jl")

# Exports
export mask_unw_igrams
export unwrapping_error
export common_scene_stack
export sbas_list
export sbas
export stack_igrams
export unw2def

# Functions
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
Deformation time series. If the network of interferograms has n scenes, n-1 maps
of the deformation (in cm) with respect to the first acquistion will be 
generated. 
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

"""
    unw2def(
        wid::Int,
        len::Int,
        ref_col::Int,
        ref_row::Int,
        intlist::String,
        max_dt::Int,
        max_err::Real,
        max_atm::Real;
        path_mask::String="",
        fraction_ram::Real=0.1,
        wavelength::Real=5.55
    )

Converts a network of unwrapped interferograms into deformation time series. 
This is an orchestrator for the following routines:
    1. apply mask and calibrate unwrapped interferograms
    2. compute unwrapping error for each interferogram
    3. estimate atmospheric delays for each scene using common scene stacking
    4. impose temporal baseline, unwrapping error, and atm. delay thresholds
    5. run SBAS and stacking on the resulting subset of data

# Positional Arguments
- `wid`: `Int` width of the unwrapped interferograms (in pixels)
- `len`: `Int` lenght of the unwrapped interferograms (in pixels)
- `ref_col`: `Int` column index of reference location
- `ref_row`: `Int` row index of reference location
- `intlist`: `String` path to file containing the list of interferograms
- `max_dt`: `Int` maximum temporal baseline (in days) 
- `max_err`: `Real` maximum total unwrapping error (in rad)
- `max_atm`: `Real` maximum atmospheric noise (in rad)

# Optional Arguments
- `path_mask`: `String` path to pixel mask (default is "" for no mask)
- `fraction_ram`: `Real` fraction of RAM to allocate to process (default is 0.1)
- `wavelength`: `Real` sensor's wavelength (in cm, default is 0.5)

# Returns
* LOS deformation time series [cm]: if the network of interferograms has n scenes, 
n-1  maps of the deformation with respect to the first acquistion will 
be generated.
* LOS deformation rate [cm/yr]: 
    1. linear SBAS solution (`sbas_rate.un`)
    2. from stacking (`stack_rate.un`)
"""
function unw2def(
    wid::Int,
    len::Int,
    ref_col::Int,
    ref_row::Int,
    intlist::String,
    max_dt::Int,
    max_err::Real,
    max_atm::Real;
    path_mask::String="",
    fraction_ram::Real=0.1,
    wavelength::Real=5.55
)
    # Apply mask, only if passed
    masked_list = intlist
    if path_mask != ""
        masked_list = mask_unw_igrams(wid, len, ref_col, ref_row, path_mask, intlist)
    end

    # Estimate unwrapping error and atmospheric delays
    unwerr_list = unwrapping_error(wid, len, masked_list)
    atmerr_list = common_scene_stack(wid, len, masked_list)

    # Apply constraints and run SBAS and stacking
    sbas_pairs = sbas_list(max_dt, max_err, max_atm, unwerr_list, atmerr_list)
    sbas(
        wid,
        len,
        ref_row,
        ref_col,
        sbas_pairs,
        fraction_ram=Float64(fraction_ram),
        wavelength=Float64(wavelength)
    )
    stack_igrams(
        wid,
        len,
        ref_col,
        ref_row,
        sbas_pairs,
        wavelength=Float64(wavelength)
    )
end

end

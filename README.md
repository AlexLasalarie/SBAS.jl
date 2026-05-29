# SBAS.jl

## Overview
SBAS.jl is a Julia package for forming InSAR deformation time series from
unwrapped interferograms.

### Input Data Format
This code expects unwrapped interferograms in the standard **ISCE/ENVI 2-Band
Flat Binary format**. 

The data must be a raw binary stream of 32-bit floating-point numbers(`Float32`) 
structured with **Band Interleaved** organization containing exactly two bands:
* **Band 1:** Amplitude (or Correlation)
* **Band 2:** Unwrapped Phase (in radians)

## Get started

### Clone the repository
Navigate to the desired directory and run:
```bash
git clone https://github.com/AlexLasalarie/SBAS.jl.git <your-dirname>
```

### Install dependencies
Navigate to the newly created directory and run:
```bash
julia install_sbas.jl
```

### Quick launch
To start a development session with all available CPU threads and auto-load the environment:
```bash
julia -t auto -i dev_startup.jl
```
This will load development tools such as Revise, BenchmarkTools, and JLD2, if installed in your 
global environment.

### Manual launch
Start a Julia session with the desired number of threads, `n`:
```bash
julia -t n
```
To lauch a session with all available threads, run:
```bash
julia -t auto
```
Enter Package mode by typing `]`:
```bash
pkg> activate .
```
Exit package mode by pressing Backspace and run:
```bash
julia> using SBAS
```

## Form time series

### Default run
Navigate to the directory containing the stack of unwrapped interferograms:
```bash
julia> cd("path/to/data")
```

To form the time series, run:
```bash
julia> sbas(wid, len, ref_row, ref_col, intlist)
```
Where `wid` is the width and `len` is the length of the unwrapped
interferograms, in pixels, `ref_row` and `ref_col` are the row and column
location of the reference pixel, and `intlist` is the path to the file containing
the list of unwrapped interferograms. This will create `n-1` (`n` is the number
of SAR scenes) maps of the LOS deformation with respect to the first acquisition
in cm).

### Custom parameters
Individual parameters can be tweaked.
```bash
sbas(
    wid::Int,
    len::Int,
    ref_row::Int,
    ref_col::Int,
    intlist::String;
    fraction_ram::Float64=0.1,
    wavelength::Float64=5.55
)
```

Get quick information from the Help mode:
```bash
julia> ?
help?> sbas
```

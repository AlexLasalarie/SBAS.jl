# Activate the current project environment
using Pkg
Pkg.activate(".")
if !isfile("Manifest.toml")
    @info "No Manifest found, instantiating..."
    Pkg.instantiate()
end

# ----- Load non-base pkgs (must be in dev. global path)

# Load Revise
try
    using Revise
    println("Revise loaded")
catch e
    @warn "Revise not found. Add to your global environment."
end

# Load BenchmarkTools
try
    using BenchmarkTools
    println("BenchmarkTools loaded")
catch e
    @warn "BenchmarkTools not found. Add to your global environment."
end

# Load JLD2 (for making test templates)
try
    using JLD2
    println("JLD2 loaded")
catch e
    @warn "JLD2 not found. Add to your global environment."
end

# Load Infiltrator
try
    using Infiltrator
    println("Infiltrator loaded")
catch e
    @warn "Infiltrator not found. Add to your global environment."
end

# ----- Load the main "Engine"
using SBAS

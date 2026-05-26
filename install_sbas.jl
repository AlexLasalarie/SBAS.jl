using Pkg

# Activate the project in the current directory
Pkg.activate(@__DIR__)

# Install all dependencies from Project.toml / Manifest.toml
Pkg.instantiate()

println("SBAS dependencies installed successfully.")

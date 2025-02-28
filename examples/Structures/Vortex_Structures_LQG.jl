# Example script for constructing an LQG vortex as a strcture

# Vortex structures allow us to construct a vortex solution where all
# fields, parameters and diagnostics are stored within a single structure

using QGDipoles

# Start by creating a grid

grid = CreateGrid(;
    Nx = 512,
    Ny = 512,
    Lx = [-5, 5],
    Ly = [-5, 5],
    cuda = false,
)

# Now we create our vortex solution, all possible keywords are included below

vortex = DefLQGVortex(grid;
    U = 1,
    ℓ = 1,
    R = Inf,
    β = 0,
    ActiveLayers = 1,
    H = 1,
    x₀ = [0, 0],
    α = 0,
    M = 8,
    tol = 1e-6,
    K₀ = nothing,
    a₀ = nothing,
    UseAnalytic = false,
    CalcVelocity = false,
    CalcVorticity = false,
    CalcEnergy = false,
    CalcEnstrophy = false,
)

# This vortex structure is immutable (i.e. can't be changed once created).
# However we can create a new vortex with U = 2 by modifying our previous one
# Note that the vortex solution is fully recalculated to ensure all parts of the
# solution are correct for the new parameters. Using low-level functionality is 
# required for making changes without recalculating the full solution.

vortex2 = UpdateVortex(grid, vortex; U = 2)

# We can plot these using `Plots.jl`

# using Plots
# plot(heatmap(grid, vortex.ψ), heatmap(grid, vortex2.ψ))

# Alternatively, we can define a parameter set using

params = DefLQGParams(;
    U = 1,
    ℓ = 1,
    R = 1,
    β = 1,
    ActiveLayers = 1,
    H = 1,
    x₀ = [0, 0],
    α = 0,
    M = 8,
    tol = 1e-6,
    K₀ = nothing,
    a₀ = nothing,
    UseAnalytic = false,
    CalcVelocity = false,
    CalcVorticity = false,
    CalcEnergy = false,
    CalcEnstrophy = false,
)

# We now create our vortex using 

vortex3 = DefLQGVortex(grid, params)

# We can modify our parameters and then make a new vortex

params2 = UpdateParams(params; α = π/4)
vortex4 = DefLQGVortex(grid, params2)

# We can plot these using `Plots.jl`

# plot(heatmap(grid, vortex3.ψ), heatmap(grid, vortex4.ψ))

# We can also move our vortex to a different grid

grid2 = CreateGrid(;
    Nx = 256,
    Ny = 256,
    Lx = [-5, 5],
    Ly = [-5, 5],
    cuda = false,
)

vortex5 = UpdateVortex(grid2, vortex4)

# And plotting gives a low resolution version of vortex4

# heatmap(grid2, vortex5.ψ)

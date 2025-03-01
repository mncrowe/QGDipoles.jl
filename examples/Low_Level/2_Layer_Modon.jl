# Example script for creating a 2-Layer QG modon

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1# vortex speed and radius
R = [1, 1]# Rossby radius in each layer
β = [0, 1]# background PV gradient in each layer

M = 8# number of coefficients in Zernike expansion
tol = 1e-8# maximum error in solution evaluation
cuda = false# use CuArrays for grid

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# Build and solve linear system for coefficients

λ = ℓ ./ R
μ = β * ℓ^2 / U

A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
K, a = SolveInhomEVP(A, B, c, d; tol)

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, q = Calc_ψq(a, U, ℓ, R, β, grid)

# Plot streamfunction ψ in both layers, if we have `Plots.jl` added

# using Plots
# plot(heatmap(grid, ψ, layer = 1), heatmap(grid, ψ, layer = 2))

# Example script for creating a 1-Layer QG modon

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1# vortex speed and radius
R = 1# Rossby radius in each layer
β = 1# background PV gradient in each layer

M = 8# number of coefficients in Zernike expansion
tol = 1e-8# maximum error in solution evaluation
cuda = false# use CuArrays for grid
method = :eigensolve# :eigensolve or :nlsolve

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# Build and solve linear system for coefficients

λ = ℓ / R
μ = β * ℓ^2 / U

A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, method)

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, q = Calc_ψq(grid, a; U, ℓ, R, β)
u, v = Calc_uv(grid, ψ)

# We can plot streamfunction ψ, if we have `Plots.jl` added

# using Plots
# heatmap(grid, ψ)

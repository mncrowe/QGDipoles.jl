# Example script for creating an N-Layer QG modon, currently set up for 5 layers

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1# vortex speed and radius
R = [1, 1, 2, 1, 2]# Rossby radius in each layer
β = [1, 0, 1, 2, 0]# background PV gradient in each layer
ActiveLayers = [1, 1, 0, 1, 1]# 1 => layer contains vortex region
K₀ = [5, 5, 5, 7]# initial guess for K in active layers

M = 8# number of coefficients in Zernike expansion
tol = 1e-8# maximum error in solution evaluation
cuda = false# use CuArrays for grid

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 20, 20

# Build and solve linear system for coefficients

λ = ℓ ./ R
μ = β * ℓ^2 / U

A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
A, B, c, d = ApplyPassiveLayers(A, B, c, d, ActiveLayers)

K, a = SolveInhomEVP(A, B, c, d; K₀, tol)
K, a = IncludePassiveLayers(K, a, ActiveLayers)

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, q = Calc_ψq(a, U, ℓ, R, β, grid)
u, v = Calc_uv(ψ, grid)

# Plot streamfunction ψ in layer 1, if we have `Plots.jl` added

# using Plots
# heatmap(grid, ψ, layer = 1)

# Example script for creating a 1-Layer QG modon

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1# vortex speed and radius
R = 1# Rossby radius in each layer
β = 1# background PV gradient in each layer

M = 8# number of coefficients in Zernike expansion
tol = 1e-8# maximum error in solution evaluation
cuda = false# use CuArrays for grid
method = 0# 0; eigensolve/nlsolve, 1; nlsolve

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# Build and solve linear system for coefficients

λ = ℓ / R
μ = β * ℓ^2 / U

A, B, c, d = BuildLinSys(M, λ, μ; tol)
K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, method)

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, q = Calc_ψq(a, U, ℓ, R, β, grid)
u, v = Calc_uv(ψ, grid)

# Plot streamfunction ψ

using Plots

heatmap(
    grid.x,
    grid.y,
    transpose(ψ[:, :, 1]);
    colormap = :balance,
    aspect_ratio = 1,
    xlims = (-Lx / 2, Lx / 2),
    ylims = (-Ly / 2, Ly / 2),
    xlabel = "x",
    ylabel = "y",
)

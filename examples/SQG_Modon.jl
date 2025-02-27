# Example script for creating an SQG QG modon
# Note: SQG modons generally require more coefficient for convergence
# 	hence a larger M. This is compensated by the faster calculation
#	of the terms in the SQG linear system.

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1# vortex speed and radius
R = [Inf, Inf]# Baroclinic and Barotropic Rossby radii
β = 0# background PV gradient in the interior

M = 20# number of coefficients in Zernike expansion
tol = 1e-6# maximum error in solution evaluation
cuda = false# use CuArrays for grid
method = 0# 0; eigensolve/nlsolve, 1; nlsolve
sqg = true# functions use SQG functionality

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# Build and solve linear system for coefficients

λ = ℓ ./ R
μ = β * ℓ^2 / U

A, B, c, d = BuildLinSys(M, λ, μ; tol, sqg)
K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, method, sqg)

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, b = Calc_ψb(a, U, ℓ, R, β, grid)
u, v = Calc_uv(ψ, grid)

# Plot surface bouyancy b

using Plots

heatmap(
    grid.x,
    grid.y,
    transpose(b);
    colormap = :balance,
    aspect_ratio = 1,
    xlims = (-Lx / 2, Lx / 2),
    ylims = (-Ly / 2, Ly / 2),
    xlabel = "x",
    ylabel = "y",
)

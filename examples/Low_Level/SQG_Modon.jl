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
method = :eigensolve# :eigensolve or :nlsolve
m = 1# exponent of K in linear system, 1 for SQG

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# Build and solve linear system for coefficients

λ = ℓ ./ R
μ = β * ℓ^2 / U

A, B, c, d = BuildLinSysSQG(M, λ, μ; tol)
K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, method, m)

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, b = Calc_ψb(grid, a; U, ℓ, R, β)
u, v = Calc_uv(grid, ψ)

# Plot surface bouyancy b, if we have `Plots.jl` added

# using Plots
# heatmap(grid, b)

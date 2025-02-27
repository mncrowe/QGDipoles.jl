# Example script for creating a 3-Layer QG modon

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1# vortex speed and radius
R = [1, 1, 1]# Rossby radius in each layer
β = [0, 0, 1]# background PV gradient in each layer
ActiveLayers = [0, 1, 0]# 1 => layer contains vortex region
x₀ = [5, 5]# location of vortex center

M = 8# number of coefficients in Zernike expansion
tol = 1e-8# maximum error in solution evaluation
cuda = false# use CuArrays for grid

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = [0, 10], [0, 10]

# Build and solve linear system for coefficients

λ = ℓ ./ R
μ = β * ℓ^2 / U

A, B, c, d = BuildLinSys(M, λ, μ; tol)
A, B, c, d = ApplyPassiveLayers(A, B, c, d, ActiveLayers)

K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol)
K, a = IncludePassiveLayers(K, a, ActiveLayers)

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, q = Calc_ψq(a, U, ℓ, R, β, grid, x₀)

# Plot streamfunction ψ in layer 2

using Plots

heatmap(
    grid.x,
    grid.y,
    transpose(ψ[:, :, 2]);
    colormap = :balance,
    aspect_ratio = 1,
    xlims = Lx,
    ylims = Ly,
    xlabel = "x",
    ylabel = "y",
)

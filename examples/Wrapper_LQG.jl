# This example demonstates how to use the wrapper function `CreateModonLQG`
# This script is currently set up to match `2_Layer_Modon.jl`

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1			# vortex speed and radius
R = [1, 1]			# Rossby radius in each layer
β = [0, 1]			# background PV gradient in each layer
ActiveLayers = [1, 1]		# 1 => layer contains vortex region
x₀ = [0, 0]			# position of vortex center

M = 8				# number of coefficients in Zernike expansion
tol = 1e-8			# maximum error in solution evaluation
cuda = false			# use CuArrays for grid
K₀, a₀ = [4, 4], Nothing	# guesses for K and a

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# create grid

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

# create modon solution

ψ, q, K, a = CreateModonLQG(grid, M, U, ℓ, R, β, ActiveLayers, x₀; K₀, a₀, tol)
# ψ, q = CreateModonLQG(grid, M, U, ℓ, R, β, ActiveLayers, x₀; K₀, a₀, tol) # fields only

# Plot streamfunction ψ in layer 1

using Plots

heatmap(grid.x, grid.y, transpose(ψ[:, :, 1]);
	colormap = :balance,
	aspect_ratio = 1,
	xlims = (-Lx/2, Lx/2),
	ylims = (-Ly/2, Ly/2),
	xlabel = "x",
	ylabel = "y")
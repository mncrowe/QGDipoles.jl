# This example demonstates how to use the wrapper function `CreateModonSQG`
# This script is currently set up to match `SQG_Modon.jl`

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1			# vortex speed and radius
R = [Inf, Inf]			# Baroclinic and Barotropic Rossby radii
β = 0				# background PV gradient in the interior
x₀ = [0, 0]			# position of vortex center

M = 20				# number of coefficients in Zernike expansion
tol = 1e-6			# maximum error in solution evaluation
cuda = false			# use CuArrays for grid
K₀, a₀ = 4, Nothing		# guesses for K and a

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# create grid

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

# create modon solution

ψ, b, K, a = CreateModonSQG(grid, M, U, ℓ, R, β, x₀; K₀, a₀, tol)
# ψ, b = CreateModonSQG(grid, M, U, ℓ, R, β, x₀; K₀, a₀, tol) # fields only

Nothing
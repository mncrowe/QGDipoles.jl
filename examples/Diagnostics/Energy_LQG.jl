# Example script for calculating the energy and enstrophy for 2-Layer QG modon

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

# Create modon solution and grid

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, q, K, a = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

# Calculate kinetic and potential energy, we've used layer depth 
# H = R (assuming g'/f^2 = 1 in nondimensional units)

KE, PE = EnergyLQG(grid, ψ; R, H = R .^ 2)

# Calculate the enstrophy

Q = EnstrophyLQG(grid, q; H = R .^ 2)

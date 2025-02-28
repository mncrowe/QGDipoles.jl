# Example script for calculating the energy of an SQG modon

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1# vortex speed and radius
R = [Inf, Inf]# Baroclinic and Barotropic Rossby radii
β = 0# background PV gradient in the interior

M = 20# number of coefficients in Zernike expansion
tol = 1e-6# maximum error in solution evaluation
cuda = false# use CuArrays for grid

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# create modon solution and grid

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, b, _, _ = CreateModonSQG(grid, M, U, ℓ, R, β; tol)

# calculate the domain integrated energy and surface potential energy

E, SPE = EnergySQG(grid, ψ, b, R[2])

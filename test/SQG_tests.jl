"""
Tests for the SQG functions from QGDipoles.jl

"""

"""
Function: `TestSQG_K()`

Tests the value of K for an SQG vortex against known values.

"""
function TestSQG_K()

    # Define parameters

    M = 13
    tol = 1e-6
    K₀ = 4
    method = 1
    sqg = true

    U, ℓ = 1, 1
    R = [Inf, Inf]
    β = 0

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    # Solve linear system

    A, B, c, d = BuildLinSys(M, λ, μ; tol, sqg)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method, sqg)

    return abs(K[1] - 4.12126) < 1e-5

end

"""
Function: `TestSQG_v(cuda)`

Tests the value of maximum(v) for an SQG vortex against known values.
This function tests the velocity calculation on the CPU and GPU grids.

"""
function TestSQG_v(cuda)

    # Define parameters

    M = 13
    tol = 1e-6
    K₀ = 4
    method = 1
    sqg = true

    U, ℓ = 1, 1
    R = [Inf, Inf]
    β = 0

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    Nx, Ny = 128, 128
    Lx, Ly = 10, 10

    # Solve linear system

    A, B, c, d = BuildLinSys(M, λ, μ; tol, sqg)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method, sqg)

    # Create solution

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
    ψ, b = Calc_ψb(a, U, ℓ, R, β, grid)

    # Calculate velocities	

    u, v = Calc_uv(ψ, grid)

    return abs(maximum(v) - 3.54912) < 1e-5

end

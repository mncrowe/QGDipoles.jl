"""
Tests for the wrapper functions from QGDipoles.jl. These tests include tests for
`CreateLCD` and `CreateLRD`.

"""

"""
Function: `TestWrapperLQG(U, ℓ, R, β; cuda)`

Tests the LQG wrapper by comparing against the results without a wrapper.

"""
function TestWrapperLQG(
    U::Number,
    ℓ::Number,
    R::Union{Number,Vector},
    β::Union{Number,Vector};
    cuda::Bool,
)

    # Set numerical parameters

    M = 8
    tol = 1e-8

    # Create grid

    Nx, Ny = 256, 256
    Lx, Ly = 20, 20

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

    # Calculate solution without wrapper

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; tol)
    ψ₁, _ = Calc_ψq(grid, a; U, ℓ, R, β)

    # Calculate solution with wrapper

    ψ₂, _, _, _ = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

    return maximum(abs.(ψ₁ - ψ₂)) < 1e-10

end

"""
Function: `TestWrapperSQG(U, ℓ, R, β; cuda)`

Tests the SQG wrapper by comparing against the results without a wrapper.

"""
function TestWrapperSQG(U::Number, ℓ::Number, R::Vector, β::Number; cuda::Bool)

    # Set numerical parameters

    M = 8
    tol = 1e-8

    # Create grid

    Nx, Ny = 256, 256
    Lx, Ly = 20, 20

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

    # Calculate solution without wrapper

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    A, B, c, d = BuildLinSysSQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; tol, m = 1)
    ψ₁, _ = Calc_ψb(grid, a; U, ℓ, R, β)

    # Calculate solution with wrapper

    ψ₂, _, _, _ = CreateModonSQG(grid; U, ℓ, R, β, M, tol)

    return maximum(abs.(ψ₁ - ψ₂)) < 1e-10

end

"""
Function: `TestLCD(U, ℓ; cuda)`

Tests the (analytic) `CreateLCD` solution by comparing against the numerical results.
Note: The LCD decays very slowly so these results actually differ in the far field as
the numerical method enforces periodicity through the Fourier transforms whereas the
analytical result decays as 1/r so remains O(0.1) at y ~ 10.

"""
function TestLCD(U::Number, ℓ::Number; cuda::Bool)

    # Set numerical parameters

    M = 8
    tol = 1e-8

    # Create grid

    Nx, Ny = 256, 256
    Lx, Ly = 20, 20

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

    # Create LCD using analytical result

    ψ₁, _, _ = CreateLCD(grid; U, ℓ)

    # Create LCD using numerical method

    ψ₂, _, _, _ = CreateModonLQG(grid; U, ℓ, R = Inf, β = 0, M, tol)

    return maximum(abs.(ψ₁ - ψ₂)) < U * ℓ * 1.5e-1

end

"""
Function: `TestLRD(U, ℓ, R, β; cuda)`

Tests the (semi-analytic) `CreateLRD` solution by comparing against the numerical results.
Note: the LRD generally decays much faster than the LCD as it is O(K_1(p*r)) in the far field.
Therefore we expect less difference between the two methods compared with the LCD results.

"""
function TestLRD(U::Number, ℓ::Number, R::Number, β::Number; cuda::Bool)

    # Set numerical parameters

    M = 8
    tol = 1e-8

    # Create grid

    Nx, Ny = 256, 256
    Lx, Ly = 20, 20

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

    # Create LRD using semi-analytical result

    ψ₁, _, _ = CreateLRD(grid; U, ℓ, R, β)

    # Create LRD using numerical method

    ψ₂, _, _, _ = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

    return maximum(abs.(ψ₁ - ψ₂)) < U * ℓ * 2e-3

end

"""
Tests for the SQG functions from QGDipoles.jl

"""

"""
    TestSQG_K()

Tests the value of K for an SQG vortex against known values.

"""
function TestSQG_K()

    # Define parameters

    M = 13
    tol = 1e-6
    K₀ = 4
    method = 1

    U, ℓ = 1, 1
    R = [Inf, Inf]
    β = 0

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    # Solve linear system

    A, B, c, d = BuildLinSysSQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method, m = 1)

    return abs(K[1] - 4.12126) < 1e-5

end

"""
    TestSQG_v(grid)

Tests the value of maximum(v) for an SQG vortex against known values.
This function tests the velocity calculation on the CPU and GPU grids.

"""
function TestSQG_v(grid)

    # Define parameters

    M = 13
    tol = 1e-6
    K₀ = 4
    method = 1

    U, ℓ = 1, 1
    R = [Inf, Inf]
    β = 0

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    # Solve linear system

    A, B, c, d = BuildLinSysSQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method, m = 1)

    # Create solution

    ψ, b = Calc_ψb(grid, a; U, ℓ, R, β)

    # Calculate velocities	

    u, v = Calc_uv(grid, ψ)

    return abs(maximum(v) - 3.54912) < 1e-5

end

"""
    TestWrapperSQG(grid; U, ℓ, R, β)

Tests the SQG wrapper by comparing against the results without a wrapper.

"""
function TestWrapperSQG(grid; U::Number, ℓ::Number, R::Vector, β::Number)

    # Set numerical parameters

    M = 8
    tol = 1e-8

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
    TestSQGVortex(grid)

Tests the construction of SQG vortex parameters and structures by building a vortex with
non-default parameters

"""
function TestSQGVortex(grid)

    vortex = DefSQGVortex(
        grid;
        U = 0.95,
        ℓ = 1.05,
        R = [1.1, 0.9],
        β = 0.8,
        x₀ = [0.01, -0.01],
        α = 0.1,
        M = 10,
        tol = 1e-5,
        K₀ = 4.05,
        a₀ = nothing,
        CalcVelocity = true,
        CalcVorticity = true,
        CalcEnergy = true,
    )

    params = UpdateParams(vortex.params; U = 1.1)

    vortex2 = UpdateVortex(grid, vortex; ℓ = 1)

    return true

end

"""
    TestSQGEnergy(grid)

Tests the values of E and SPE for an SQG model against known values.

"""
function TestSQGEnergy(grid)

    # Set parameters

    U, ℓ, R, β = 1, 1, [Inf, Inf], 0
    M = 8
    tol = 1e-8

    # Create LRD using numerical method

    ψ, b, _, _ = CreateModonSQG(grid; U, ℓ, R, β, M, tol)

    E, SPE = EnergySQG(grid, ψ, b; R′ = R[2])

    TestE = E[1] - 9.750526656500448 < 1e-6
    TestSPE = SPE[1] - 30.137121991640193 < 1e-6

    return TestE & TestSPE

end

"""
Tests for the LQG functiond for QGDipoles.jl

"""

"""
    Test1QG_K(U, ℓ, R, β)

Tests the value of K for a 1-layer QG solution by comparing with the result obtained by solving
the root finding problem of Larichev & Reznik.

"""
function Test1QG_K(U::Number, ℓ::Number, R::Number, β::Number)

    # Define parameters

    M = 8
    tol = 1e-6
    K₀ = 4
    method = 1

    λ = ℓ / R
    μ = β * ℓ^2 / U

    # Create linear system

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K₁, a = SolveInhomEVP(A, B, c, d; K₀, tol, method)

    # Calculate semi-analytic result using root finding

    β′ = β + U / R^2
    p = sqrt(β′ / U)

    if p == 0

        K₂ = 3.83170597020751231561443589 / ℓ

    else

        # Define Bessel functions and derivatives

        J1(x) = besselj(1, x)
        J1p(x) = (besselj(0, x) - besselj(2, x)) / 2
        K1(x) = besselk(1, x)
        K1p(x) = (-besselk(0, x) - besselk(2, x)) / 2

        # Define a function f(x), K is related to the zeros of f

        f(x) = @. x * J1p(x) - (1 + x^2 / (p^2 * ℓ^2)) * J1(x) +
           x^2 * J1(x) * K1p(p * ℓ) / (p * ℓ * K1(p * ℓ))

        # Solve f(x) = 0 and calculate K

        K′ = nlsolve(f, [3.83170597020751231561443589]).zero[1] / ℓ
        K₂ = ℓ * sqrt(K′^2 + 1 / R^2)

    end

    return abs(K₁[1] - K₂) < 1e-5

end

"""
    Test1QG_ψ(cuda)

Tests the value of maximum value of ψ for the LCD against a known value.

"""
function Test1QG_ψ(cuda::Bool)

    # Define parameters

    M = 8
    tol = 1e-6
    K₀ = 4
    method = 1

    λ, μ = 0, 0

    Nx, Ny = 128, 128
    Lx, Ly = 10, 10

    # Build linear system

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method)

    # Calculate solution

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda = cuda)
    ψ, q = Calc_ψq(a, 1, 1, Inf, 0, grid)

    return abs(maximum(ψ) - 1.28110) < 1e-5

end

"""
    Test2QG_K()

Tests the value of (K_1, K_2) for a 2-layer vortex against known values.

"""
function Test2QG_K()

    # Define parameters

    M = 8
    tol = 1e-6
    K₀ = [4, 4]
    method = 1

    U, ℓ = 1, 1
    R = [1, 1]
    β = [0, 1]

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    # Build linear system

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K₁, a = SolveInhomEVP(A, B, c, d; K₀, tol, method)

    # Set known value for K

    K₂ = [3.8002227321789492, 3.949863664432825]

    return maximum(abs.(@. K₁ - K₂')) < 1e-5

end

"""
    Test2QG_PVinv(cuda)

Tests if the 3D PV inversion array has been built correctly to avoid Inf/NaN values in ψ.
Common error in (0, 0) wavenumber due to divide by 0 issues. Due to differences in function
support, these arrays are constructed differently on the CPU and GPU.

"""
function Test2QG_PVinv(cuda)

    # Define parameters

    M = 8
    tol = 1e-6
    K₀ = [4, 4]
    method = 1

    U, ℓ = 1, 1
    R = [Inf, 1]
    β = [0, 1]

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    Nx, Ny = 128, 128
    Lx, Ly = 10, 10

    # Build linear system

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method)

    # Calculate solution

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
    ψ, q = Calc_ψq(a, U, ℓ, R, β, grid)

    # Look for NaN and Inf values

    No_NaNs_ψ = ~maximum(isnan.(ψ))
    No_Infs_ψ = ~maximum(isinf.(ψ))
    No_NaNs_q = ~maximum(isnan.(q))
    No_Infs_q = ~maximum(isinf.(q))

    return (No_NaNs_ψ & No_Infs_ψ) & (No_NaNs_q & No_Infs_q)

end

"""
    TestLQG_ActiveLayers()

Tests if active and passive layers are being correctly applied. This is done by building a three
layer system and comparing the results obtained by removing the passive layers with the results
obtained by looking for negative eigenvalues for passive layers (K^2 = -β/U). We also check if
the reduced linear system has the correct size, this should be M x A, where A is number of active
layers, here A = 1. 

"""
function TestLQG_ActiveLayers()

    # Set parameters

    M = 8
    tol = 1e-8
    K₀ = [0, 4, 0]# use 0 as guess for K^2 = -β/U values as close enough to converge
    warn = false# suppress warning that solution includes passive layers

    U, ℓ = 1, 1
    R = [1, 1, 1]
    β = [0, 0, 1]
    ActiveLayers = [0, 1, 0]

    # Build and solve full linear system

    λ = ℓ ./ R
    μ = β * ℓ^2 / U

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K₁, _ = SolveInhomEVP(A, B, c, d; K₀, tol, warn)

    # Build and solve reduced linear system

    A, B, c, d = ApplyPassiveLayers(A, B, c, d, ActiveLayers)
    K₂, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol)
    K₂, _ = IncludePassiveLayers(K₂, a, ActiveLayers)

    # Check K values match

    If_K = maximum(abs.(K₁ - K₂)) < 1e-6

    # Check system size

    If_N = size(A) == (M, M)

    return If_K & If_N

end

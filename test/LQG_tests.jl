"""
Tests for the LQG functions for QGDipoles.jl

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
function Test1QG_ψ(grid)

    # Define parameters

    M = 8
    tol = 1e-6
    K₀ = 4
    method = 1

    λ, μ = 0, 0

    # Build linear system

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method)

    # Calculate solution

    ψ, q = Calc_ψq(grid, a; U = 1, ℓ = 1, R = Inf, β = 0)

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
    Test2QG_PVinv(grid)

Tests if the 3D PV inversion array has been built correctly to avoid Inf/NaN values in ψ.
Common error in (0, 0) wavenumber due to divide by 0 issues. Due to differences in function
support, these arrays are constructed differently on the CPU and GPU.

Also test length(β) < N case through use of β = 0 for N = 2.

"""
function Test2QG_PVinv(grid)

    # Define parameters

    M = 8
    tol = 1e-6
    K₀ = [4, 4]
    method = 1

    U, ℓ = 1, 1
    R = [Inf, 1]
    β = 0

    λ = ℓ ./ R
    μ = [β, β] * ℓ^2 / U

    # Build linear system

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    K, a = SolveInhomEVP(A, B, c, d; K₀, tol, method)

    # Calculate solution

    ψ, q = Calc_ψq(grid, a; U, ℓ, R, β)

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

"""
    TestWrapperLQG(grid; U, ℓ, R, β)

Tests the LQG wrapper by comparing against the results without a wrapper.

"""
function TestWrapperLQG(
    grid;
    U::Number,
    ℓ::Number,
    R::Union{Number,Vector},
    β::Union{Number,Vector},
)

    # Set numerical parameters

    M = 8
    tol = 1e-8

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
    TestLCD(U, ℓ; cuda)

Tests the (analytic) `CreateLCD` solution by comparing against the numerical results.
Note: The LCD decays very slowly so these results actually differ in the far field as
the numerical method enforces periodicity through the Fourier transforms whereas the
analytical result decays as 1/r so remains O(0.1) at y ~ 10.

"""
function TestLCD(
    grid;
    U::Number,
    ℓ::Number,
)

    # Set numerical parameters

    M = 8
    tol = 1e-8

    # Create LCD using analytical result

    ψ₁, _, _ = CreateLCD(grid; U, ℓ)

    # Create LCD using numerical method

    ψ₂, _, _, _ = CreateModonLQG(grid; U, ℓ, R = Inf, β = 0, M, tol)

    return maximum(abs.(ψ₁ - ψ₂)) < U * ℓ * 1.5e-1

end

"""
    TestLRD(U, ℓ, R, β; cuda)

Tests the (semi-analytic) `CreateLRD` solution by comparing against the numerical results.
Note: the LRD generally decays much faster than the LCD as it is O(K_1(p*r)) in the far field.
Therefore we expect less difference between the two methods compared with the LCD results.

"""
function TestLRD(
    grid;
    U::Number,
    ℓ::Number,
    R::Number,
    β::Number,
)

    # Set numerical parameters

    M = 8
    tol = 1e-8

    # Create LRD using semi-analytical result

    ψ₁, _, _ = CreateLRD(grid; U, ℓ, R, β)

    # Create LRD using numerical method

    ψ₂, _, _, _ = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

    return maximum(abs.(ψ₁ - ψ₂)) < U * ℓ * 2e-3

end

"""
    TestLQGVortex(grid)

Tests the construction of LQG vortex parameters and structures by building a vortex with
non-default parameters

"""
function TestLQGVortex(grid)

    vortex = DefLQGVortex(
        grid;
        U = 0.95,
        ℓ = 1.05,
        R = [1.1, 0.9],
        β = [0.8, 1.2],
        ActiveLayers = [1, 1],
        H = [0.25, 0.75],
        x₀ = [0.01, -0.01],
        α = 0.1,
        M = 6,
        tol = 1e-5,
        K₀ = [4.1, 4.2],
        a₀ = nothing,
        UseAnalytic = false,
        CalcVelocity = true,
        CalcVorticity = true,
        CalcEnergy = true,
        CalcEnstrophy = true,
    )

    params = UpdateParams(vortex.params; U = 1.1)

    vortex2 = UpdateVortex(grid, vortex; ℓ = 1)

    return true

end

"""
    TestLQGEnergy1Layer(cuda)

Tests the values of KE and PE for a 1-layer QG model against known values.

"""
function TestLQGEnergy1Layer(grid)

    # Set parameters

    U, ℓ, R, β = 1, 1, 1, 1
    M = 8
    tol = 1e-8

    # Create LRD using numerical method

    ψ, _, _, _ = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

    KE, PE = EnergyLQG(grid, ψ; R)

    TestKE = KE[1] - 12.653032138613069 < 1e-6
    TestPE = PE[1] - 2.2180534691913154 < 1e-6

    return TestKE & TestPE

end

"""
    TestLQGEnergy2Layer(grid)

Tests the values of KE and PE for a 2-layer QG model against known values.

"""
function TestLQGEnergy2Layer(grid)

    # Set parameters

    U, ℓ, R, β = 1, 1, [1, 1], [0, 1]
    M = 8
    tol = 1e-8

    # Create LRD using numerical method

    ψ, _, _, _ = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

    KE, PE = EnergyLQG(grid, ψ; R, H = [1, 1])

    TestKE1 = KE[1] - 3.5838305725995547 < 1e-6
    TestKE2 = KE[2] - 4.7052971811344300 < 1e-6
    TestPE = PE[1] - 0.0309253993785872 < 1e-6

    return TestKE1 & TestKE2 & TestPE

end


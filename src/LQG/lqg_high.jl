"""
# High-Level Functions

This file contains high-level functions. These are referred to as 'wrappers' as they
'wrap' the low-level functions into easy-to-use methods for calculating vortex solutions
without needing to understand the underlying method of linear algebra problem.

This file also includes functions which calculate the Lamb-Chaplygin Dipole and the
Larichev-Reznik Dipole using the analytic solutions for these cases.

"""

"""
    CreateModonLQG(grid; U=1, ℓ=1, R=Inf, β=0, ActiveLayers=1, x₀=[0, 0], α=0, M=8, tol=1e-6, K₀=nothing, a₀=nothing)

High level wrapper function for calculating ``ψ`` and ``q`` for the Layered QG model using given parameters

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - (`U`, `ℓ`): vortex speed and radius, Numbers (default: (`1`, `1`))
 - (`R`, `β`): Rossby radii and (y) PV gradients in each layer, Numbers or Vectors, (default: (`Inf`, `0`))
 - `ActiveLayers`: vector of 1s or 0s where 1 denotes an active layer, Number or Vector, (default: `[1,..,1]`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)
 - `α`: initial angle of vortex, Number (default: 0)
 - `M`: number of coefficient to solve for, Integer (default: `8`)
 - `tol`: error tolerance passed to `QuadGK` and `NLSolve` functions, Number (default: `1e-6`)
 - `K₀`, `a₀`: initial guesses for ``K`` and ``a``, Arrays or nothings (default: `nothing`)

Note: provide values of K₀ and a₀ for active layers ONLY.
"""
function CreateModonLQG(
    grid;
    U::Number = 1,
    ℓ::Number = 1,
    R::Union{Number,Vector} = Inf,
    β::Union{Number,Vector} = 0,
    ActiveLayers::Union{Number,Vector} = 1,
    x₀::Vector = [0, 0],
    α::Number = 0,
    M::Int = 8,
    tol = 1e-6,
    K₀::Union{Number,Array,Nothing} = nothing,
    a₀::Union{Array,Nothing} = nothing,
)

    # If ActiveLayers size does not match size of R, assume all layers are active

    if length(ActiveLayers) < length(R)

        ActiveLayers = ones(length(R))

    end

    # Define intermediate variables

    λ, μ = ℓ ./ R, β .* (ℓ^2 / U)

    # Build linear system

    A, B, c, d = BuildLinSysLQG(M, λ, μ; tol)
    A, B, c, d = ApplyPassiveLayers(A, B, c, d, ActiveLayers)

    # Solve linear system

    K, a = SolveInhomEVP(A, B, c, d; K₀, a₀, tol)
    K, a = IncludePassiveLayers(K, a, ActiveLayers)

    # Construct solution using computed coefficients

    ψ, q = Calc_ψq(grid, a; U, ℓ, R, β, x₀, α)

    return ψ, q, K, a

end

"""
    CreateLCD(grid; U=1, ℓ=1, x₀=[0, 0], α=0)

High level wrapper function for calculating ``ψ`` and ``q`` for the Lamb-Chaplygin dipole using given parameters

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - (`U`, `ℓ`): vortex speed and radius, Numbers (default: (`1`, `1`))
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)
 - `α`: initial angle of vortex, Number (default: `0`)

Note: This function uses the analytic solution for the LCD to calculate ``ψ`` and ``q``.
"""
function CreateLCD(
    grid;
    U::Number = 1,
    ℓ::Number = 1,
    x₀::Vector = [0, 0],
    α::Number = 0,
)

    # Define K as the first root of Bessel J_1(x)

    K = 3.83170597020751231561443589 / ℓ

    # Define Coefficients in far-field (A) and inside vortex (B)

    A = -U * ℓ^2
    B = 4 * U / (K * (besselj(0, K * ℓ) - besselj(2, K * ℓ)))

    # Create Cartesian and polar grids

    x, y = CartesianGrid(grid)
    r, θ = PolarGrid(x, y, x₀)

    # Calculate ψ and q using analytic result

    ψ = @. (A / r * (r >= ℓ) + (B * besselj(1, K * r) - U * r) * (r < ℓ)) * sin(θ - α)
    q = @. -K^2 * B * besselj(1, K * r) * (r < ℓ) * sin(θ - α)

    # Move result to GPU if `cuda=true` in grid

    if grid.Krsq isa CuArray

        ψ, q = CuArray(ψ), CuArray(q)

    end

    return ψ, q, [K;;]

end

"""
    CreateLRD(grid; U=1, ℓ=1, R=Inf, β=0, x₀=[0, 0], α=0)

High level wrapper function for calculating ``ψ`` and ``q`` for the Larichev-Reznik dipole using given parameters

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - (`U`, `ℓ`): vortex speed and radius, Numbers (default: (`1`, `1`))
 - (`R`, `β`): Rossby radii and (y) PV gradient, Numbers, (default: (`Inf`, `0`))
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)
 - `α`: initial angle of vortex, Number (default: `0`)

Note: This function uses the analytic solution for the LRD to calculate ``ψ`` and ``q``.
"""
function CreateLRD(
    grid;
    U::Number = 1,
    ℓ::Number = 1,
    R::Number = Inf,
    β::Number = 0,
    x₀::Vector = [0, 0],
    α::Number = 0,
)

    # Define effective β and external wavenumber p

    β′ = β + U / R^2
    p = sqrt(β′ / U)

    if p == 0

        # If parameters match the LCD, use `CreateLCD` instead

        return CreateLCD(grid, U, ℓ, x₀, α)

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
        K = ℓ * sqrt(K′^2 + 1 / R^2)

        # Define Coefficients in far-field (A) and inside vortex (B)

        A = -U * ℓ / K1(p * ℓ)
        B = p^2 * U * ℓ / (K′^2 * J1(K′ * ℓ))

        # Create Cartesian and polar grids

        x, y = CartesianGrid(grid)
        r, θ = PolarGrid(x, y, x₀)

        # Calculate ψ and q using analytic result

        ψ = @. (
            A * K1(p * r) * (r >= ℓ) +
            (B * J1(K′ * r) - U * (K′^2 + p^2) / K′^2 * r) * (r < ℓ)
        ) * sin(θ - α)
        q = @. β / U * ψ * (r >= ℓ) -
           (K^2 / ℓ^2 * ψ + (U * K^2 / ℓ^2 + β) * r * sin(θ - α)) * (r < ℓ)

        # Move result to GPU if `cuda=true` in grid

        if grid.Krsq isa CuArray

            ψ, q = CuArray(ψ), CuArray(q)

        end

        return ψ, q, [K;;]

    end

end
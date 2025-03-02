"""
# High-Level Functions

This file contains high-level functions. These are referred to as 'wrappers' as they
'wrap' the low-level functions into easy-to-use methods for calculating vortex solutions
without needing to understand the underlying method of linear algebra problem.

This file also contains functions which calculate ψ, b, q and w at specified depths using
the SQG (surface) solution for ψ.

"""

"""
    CreateModonSQG(grid; U=1, ℓ=1, R=[Inf, Inf], β=0, x₀=[0, 0], α=0, M=12, tol=1e-6, K₀=nothing, a₀=nothing)

High level wrapper function for calculating ``ψ`` and ``b`` for the SQG model using given parameters

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `M`: number of coefficient to solve for, Integer (default: `12`)
 - (`U`, `ℓ`): vortex speed and radius, Numbers (default: (`1`, `1`))
 - `R`: vector of ``[R, R']``, Vector (default: `[Inf, Inf]`)
 - `β`: beta-plane (y) PV gradient, Number (default: `0`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)
 - `α`: initial angle of vortex, Number (default: `0`)
 - `K₀`, `a₀`: initial guesses for ``K`` and ``a``, Arrays or nothings (default: `nothing`)
 - `tol`: error tolerance passed to `QuadGK` and `NLSolve` functions, Number (default: `1e-6`)

Note: Here R is the baroclinic Rossby radius, R = NH/f, and R' = R₀²/R where R₀ is
the barotropic Rossby radius, R₀ = √(gH)/f. For infinite depth, R' = g/(fN).
"""
function CreateModonSQG(
    grid;
    U::Number = 1,
    ℓ::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
    x₀::Vector = [0, 0],
    α::Number = 0,
    M::Int = 12,
    tol = 1e-6,
    K₀::Union{Number,Array,Nothing} = nothing,
    a₀::Union{Array,Nothing} = nothing,
)

    # Define intermediate variables

    λ, μ = ℓ ./ R, β .* (ℓ^2 / U)

    # Build linear system

    A, B, c, d = BuildLinSysSQG(M, λ, μ; tol)

    # Solve linear system

    K, a = SolveInhomEVP(A, B, c, d; K₀, a₀, tol, m = 1)

    # Construct solution using computed coefficients

    ψ, b = Calc_ψb(grid, a; U, ℓ, R, β, x₀, α)

    return ψ, b, K, a

end

"""
    Eval_ψ_SQG(grid, ψ; z=[0], U=1, R=[Inf, Inf], β=0)

Evaluates ``ψ`` at specified depths, ``z ∈ [-R, 0]``, for the SQG problem

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `ψ`: surface streamfunction, calculated using `Calc_ψb` or `CreateModonSQG`
 - `z`: vector of depths (default: `[0]`)
 - `U`: vortex speed, Number (default: `1`)
 - `R`: vector of ``[R, R']``, Vector (default: `[Inf, Inf]`)
 - `β`: beta-plane (y) PV gradient, Number (default: `0`)

Note: Here ``R`` is the baroclinic Rossby radius, ``R = NH/f``, and ``R' = R₀²/R`` where ``R₀`` is
the barotropic Rossby radius, ``R₀ = √(gH)/f``. For infinite depth, ``R' = g/(fN)``.
"""
function Eval_ψ_SQG(
    grid,
    ψ::Union{CuArray,Array};
    z::Vector = [0],
    U::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
)

    Nx = length(grid.x)
    ϵ = 1e-15

    # Reshape z so the z direction corresponds to the third dimension of the output

    z = reshape(z, 1, 1, :)

    # Move inputs to GPU if `cuda=true` in grid

    if grid.Krsq isa CuArray

        z, ψ = CuArray(z), CuArray(ψ)

    end

    # Fourier tranform ψ

    ψh = rfft(ψ)

    # Calculate exponents for analytic solution in Fourier space

    k₁ = @. (ϵ + sqrt(grid.Krsq + β / U)) * z
    k₂ = @. (ϵ + sqrt(grid.Krsq + β / U)) * R[1]

    # Calculate ψ in Fourier space, we divide through by exp(k₂) to prevent Inf values

    ψ₃h = @. ψh * (exp(k₁) + exp(-k₁ - 2 * k₂)) / (1 + exp(-2 * k₂))

    # Transform back to real space

    ψ₃ = irfft(ψ₃h, Nx, [1, 2])

    return ψ₃

end

"""
    Eval_q_SQG(grid, ψ; z=[0], U=1, R=[Inf, Inf], β=0)

Evaluates ``q`` at specified depths, ``z ∈ [-R, 0]``, for the SQG problem

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `ψ`: surface streamfunction, calculated using `Calc_ψb` or `CreateModonSQG`
 - `z`: vector of depths (default: `[0]`)
 - `U`: vortex speed, Number (default: `1`)
 - `R`: vector of ``[R, R']``, Vector (default: `[Inf, Inf]`)
 - `β`: beta-plane (y) PV gradient, Number (default: `0`)

Note: Here ``R`` is the baroclinic Rossby radius, ``R = NH/f``, and ``R' = R₀²/R`` where ``R₀`` is
the barotropic Rossby radius, ``R₀ = √(gH)/f``. For infinite depth, ``R' = g/(fN)``.
"""
function Eval_q_SQG(
    grid,
    ψ::Union{CuArray,Array};
    z::Vector = [0],
    U::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
)

    # Calculate q using the result q = (β / U) * ψ in the 3D domain

    q₃ = β / U * Eval_ψ_SQG(grid, ψ, z, U, R, β)

    return q₃

end

"""
    Eval_b_SQG(grid, ψ; z=[0], U=1, R=[Inf, Inf], β=0)

Evaluates ``b`` at specified depths, ``z ∈ [-R, 0]``, for the SQG problem

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `ψ`: surface streamfunction, calculated using `Calc_ψb` or `CreateModonSQG`
 - `z`: vector of depths (default: `[0]`)
 - `U`: vortex speed, Number (default: `1`)
 - `R`: vector of ``[R, R']``, Vector (default: `[Inf, Inf]`)
 - `β`: beta-plane (y) PV gradient, Number (default: `0`)

Note: Here ``R`` is the baroclinic Rossby radius, ``R = NH/f``, and ``R' = R₀²/R`` where ``R₀`` is
the barotropic Rossby radius, ``R₀ = √(gH)/f``. For infinite depth, ``R' = g/(fN)``.
"""
function Eval_b_SQG(
    grid,
    ψ::Union{CuArray,Array};
    z::Vector = [0],
    U::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
)

    Nx = length(grid.x)
    ϵ = 1e-15

    # Reshape z so the z direction corresponds to the third dimension of the output

    z = reshape(z, 1, 1, :)

    # Move inputs to GPU if `cuda=true` in grid

    if grid.Krsq isa CuArray

        z, ψ = CuArray(z), CuArray(ψ)

    end

    # Fourier tranform ψ

    ψh = rfft(ψ)

    # Calculate exponents for analytic solution in Fourier space

    k₁ = @. (ϵ + sqrt(grid.Krsq + β / U)) * z
    k₂ = @. (ϵ + sqrt(grid.Krsq + β / U)) * R[1]

    # Calculate ψ in Fourier space, we divide through by exp(k₂) to prevent Inf values

    b₃h =
        @. sqrt(grid.Krsq + β / U) * ψh * (exp(k₁) - exp(-k₁ - 2 * k₂)) / (1 + exp(-2 * k₂))

    # Transform back to real space

    b₃ = irfft(b₃h, Nx, [1, 2])

    return b₃

end

"""
    Eval_w_SQG(grid, ψ; z=[0], U=1, R=[Inf, Inf], β=0)

Evaluates N²w at specified depths, ``z ∈ [-R, 0]``, for the SQG problem using ``N²w = -J[ψ + Uy, b]``

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `ψ`: surface streamfunction, calculated using `Calc_ψb` or `CreateModonSQG`
 - `z`: vector of depths (default: `[0]`)
 - `U`: vortex speed, Number (default: `1`)
 - `R`: vector of ``[R, R']``, Vector (default: `[Inf, Inf]`)
 - `β`: beta-plane (y) PV gradient, Number (default: `0`)

Note: Here ``R`` is the baroclinic Rossby radius, ``R = NH/f``, and ``R' = R₀²/R`` where ``R₀`` is
the barotropic Rossby radius, ``R₀ = √(gH)/f``. For infinite depth, ``R' = g/(fN)``.

Note: this function is not accurate at the surface as ``∇b`` is discontinuous there.
Instead use ``w = -U∂η/∂x`` where ``η = fψ/g`` is the surface elevation, or ``w = 0`` if ``R' = ∞``.
"""
function Eval_w_SQG(
    grid,
    ψ::Union{CuArray,Array};
    z::Vector = [0],
    U::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
)
    # Calculate ψ and b at given depth

    ψ₃ = Eval_ψ_SQG(grid, ψ, z, U, R, β)
    b₃ = Eval_b_SQG(grid, ψ, z, U, R, β)

    # Calculate velocities and buoyancy gradients

    u₃, v₃ = Calc_uv(ψ₃, grid)
    b₃x, b₃y = Calc_∇(b₃, grid)

    # Evaluate N²w = (U-u)∂b/∂x - v∂b/∂x

    w = @. -((u₃ - U) * b₃x + v₃ * b₃y)

    return w

end

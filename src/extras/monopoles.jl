"""
While primarily intended for creating dipolar vortex solutions, this package contains
additional functions which can be used to generate some simple monopolar vortex solutions
too. These functions are included in this file.

An additional use of these functions is to build 'riders' as discussed by Kizner et. al.
2003 which are monopolar vortices superimposed on our dipole solutions.
"""

"""
    CreateRankine(grid; ℓ=1, Γ=2π, x₀=[0, 0])

Calculates the Rankine vortex for a 1 layer system. This vortex appears as a
point vortex in the far field but consists of solid body rotation within the
region ``r < ℓ``.

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `ℓ`: vortex speed and radius, Numbers (default: `1`)
 - `Γ`: vortex circulation (default: `2π`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)

Note: This function outputs ``(u, v)`` directly since the solution has discontinuous velocity at
the vortex boundary, ``r = ℓ``, so derivatives evaluated with Fourier transforms exhibit Gibbs
phenomenon.
"""
function CreateRankine(
    grid;
    ℓ::Number = 1,
    Γ::Number = 2π,
    x₀::Vector = [0, 0],
)

    # Create Cartesian and polar grids

    x, y = CartesianGrid(grid)
    r, θ = PolarGrid(x, y, x₀)

    # Calculate ψ and q using analytic result

    ψ = @. Γ / (2π) * (log(r / ℓ) * (r >= ℓ) + 1 / 2 * (r^2 / ℓ^2 - 1) * (r < ℓ))
    q = @. Γ / (2π) * (2 / ℓ^2) * (r < ℓ)
    u = @. -Γ / (2π) * (1 / r * (r >= ℓ) + r / ℓ^2 * (r < ℓ)) * sin(θ)
    v = @. Γ / (2π) * (1 / r * (r >= ℓ) + r / ℓ^2 * (r < ℓ)) * cos(θ)

    # Move result to GPU if `cuda=true` in grid

    if grid.Krsq isa CuArray

        ψ = CuArray(ψ)
        q = CuArray(q)
        u = CuArray(u)
        v = CuArray(v)

    end

    return ψ, q, u, v

end

"""
    Create1LMonopole(grid; ℓ=1, Γ=2π, R=Inf, x₀=[0, 0])

Calculates a monopolar vortex satisfying a Long's model assumption ``q = F(ψ)``
where ``q = [∇²-1/R²]ψ``. We take ``F(z) = -(K²+1/R²)(z-z₀)`` for ``r < ℓ`` and ``F(z) = 0``
for ``r > ℓ`` and ``z₀ = ψ(r=ℓ)``. These solutions exist only on an f-plane ``(β = 0)``.

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `ℓ`: vortex speed and radius, Numbers (default: `1`)
 - `Γ`: vortex circulation (default: `2π`)
 - `R`: Rossby radius (default: `Inf`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)

Note: This vortex has a continuous vorticity distribution so calculating (u, v)
from ψ with Fourier transforms will work. This function outputs (u, v) from the
analytical expressions for consistency with `CreateRankine`.
"""
function Create1LMonopole(
    grid;
    ℓ::Number = 1,
    Γ::Number = 2π,
    R::Number = Inf,
    x₀::Vector = [0, 0],
)

    # Define Bessel functions and derivatives

    J0(x) = besselj(0, x)
    J0p(x) = -besselj(1, x)
    K0(x) = besselk(0, x)
    K0p(x) = -besselk(1, x)

    if R == Inf

        # Use analytic results for coefficients

        K = 2.40482555769577 / ℓ
        B = Γ / (2π)
        A = B / (J0p(K * ℓ) * ℓ)

        # Create Cartesian and polar grids

        x, y = CartesianGrid(grid)
        r, θ = PolarGrid(x, y, x₀)

        # Calculate ψ and q using analytic result

        ψ = @. A * J0(K * r) * (r < ℓ) + B * log(r / ℓ) * (r >= ℓ)
        q = @. -K^2 * A * J0(K * r) * (r < ℓ)
        u = @. -(A * K * J0p(K * r) * (r < ℓ) + B / r * (r >= ℓ)) * sin(θ)
        v = @. (A * K * J0p(K * r) * (r < ℓ) + B / r * (r >= ℓ)) * cos(θ)

    else

        # Define a function f(x), K is given by the zeros of f

        f(x) = @. J0p(x) * K0(ℓ / R) + x * (R / ℓ) * J0(x) * K0p(ℓ / R)

        # Solve f(x) = 0 for x = Kℓ and set coefficients

        K = nlsolve(f, [2.40482555769577]).zero[1] / ℓ
        B = Γ / (2π)
        A = -B * K0(ℓ / R) / (K^2 * R^2 * J0(K * ℓ))

        # Create Cartesian and polar grids

        x, y = CartesianGrid(grid)
        r, θ = PolarGrid(x, y, x₀)

        # Calculate ψ and q using analytic result

        ψ = @. (A * J0(K * r) + B * (1 + 1 / (K^2 * R^2)) * K0(ℓ / R)) * (r < ℓ) +
           B * K0(r / R) * (r >= ℓ)
        q = @. -(K^2 + 1 / R^2) * (A * J0(K * r) + B / (K^2 * R^2) * K0(ℓ / R)) * (r < ℓ)
        u = @. -(A * K * J0p(K * r) * (r < ℓ) + B / R * K0p(r / R) * (r >= ℓ)) * sin(θ)
        v = @. (A * K * J0p(K * r) * (r < ℓ) + B / R * K0p(r / R) * (r >= ℓ)) * cos(θ)

    end

    # Move result to GPU if `cuda=true` in grid

    if grid.Krsq isa CuArray

        ψ = CuArray(ψ)
        q = CuArray(q)
        u = CuArray(u)
        v = CuArray(v)

    end

    return ψ, q, u, v

end

"""
    InvertVorticity1LQG(grid, q; R=Inf)

This function inverts the potential vorticity relation ``q = [∇²-1/R²]ψ`` for 1-layer QG

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `q`: potential vorticity field, Array
 - `R`: Rossby radius (default: `Inf`)

Note: This function is designed to be used for creating periodic streamfunctions using
the vorticity fields generated by `CreateMonopole`. It does not support multi-layer QG
and is only valid on an f-plane (β = 0).
"""
function InvertVorticity1LQG(
    grid,
    q::Union{CuArray, Array};
    R::Number = Inf,
)

    Nx = size(q, 1)

    # Fourier transform q

    qh = rfft(q)

    # If R = Inf we must zero the (0, 0) mode to prevent div by 0 errors

    if R == Inf

        CUDA.@allowscalar(qh[1, 1] = 0)

    end

    # Calculate ψ in Fourier space and transform back to real space

    ψh = @. -qh / (grid.Krsq + 1 / R^2)
    ψ = irfft(ψh, Nx)

    return ψ

end

"""
    CreateLQGMonopole(grid; ℓ=1, E=1, R=Inf, x₀=[0, 0])

Calculates a monopolar vortex in the LQG model using a numerical approach.
We assume that ``qⱼ + βⱼ = Fⱼ(ψⱼ + Uy)`` and write ``Fⱼ(z) = -Kⱼ² z + Eⱼ``. Expanding
the expression gives ``qⱼ + βⱼ = -Kⱼ²(ψⱼ + Uy) + Eⱼ`` which by linearity can be
split into a dipole equation ``qⱼ + βⱼ = -Kⱼ²(ψⱼ + Uy)`` and a monopole equation
``qⱼ = Eⱼ``. Outside the vortex, we take ``qⱼ = 0``.

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `ℓ`: vortex speed and radius, Numbers (default: `1`)
 - `E`: vector of ``Eⱼ`` values, Number or Vector (default: `[1, ... , 1]`)
 - `R`: Rossby radius (default: `Inf`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)

"""
function CreateLQGMonopole(
    grid;
    ℓ::Number = 1,
    E::Union{Vector, Number} = 1,
    R::Union{Vector, Number} = Inf,
    x₀::Vector = [0, 0],
)

    # Get coordinates

    x, y = CartesianGrid(grid)
    r, _ = PolarGrid(x, y, x₀)

    # Set PV in each layer

    q = (r .<= ℓ) .* reshape(E, 1, 1, :)

    if grid.Krsq isa CuArray

        q = CuArray(q)

    end

    # calculate ψ

    ψ = InvertVorticityLQG(grid, q, R)

    return ψ, q

end

"""
    InvertVorticityLQG(grid, q; R=Inf)

This function inverts the potential vorticity relation ``q = ΔN ψ`` for the LQG model

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `q`: potential vorticity field, Array
 - `R`: Rossby radius, Number or Vector (default: `Inf`)

"""
function InvertVorticityLQG(
    grid,
    q::Union{CuArray, Array};
    R::Union{Vector, Number} = Inf,
)

    Nx, Ny = size(q)
    N = Int(length(q) / (Nx * Ny))

    # Define ΔN operator, where ΔN ψ = q

    ΔN = ΔNCalc(grid.Krsq, R, 0)

    # Invert ΔN

    ΔN_i = stack(inv, eachslice(ΔN, dims = (3, 4)))
    CUDA.@allowscalar(ΔN_i[:, :, 1, 1] .= 0)

    # Calculate qh and define ψh

    qh = rfft(q, [1, 2])
    ψh = 0 .* qh

    # Calculate ψh from qh in Fourier space

    for n = 1:N
        for j = 1:N
            ψh[:, :, n] .+= ΔN_i[n, j, :, :] .* qh[:, :, j]
        end
    end

    ψ = irfft(ψh, Nx, [1, 2])

    return ψ

end


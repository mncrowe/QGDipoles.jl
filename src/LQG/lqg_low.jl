"""
# Low-Level Functions

This file contains functions for calculating LQG vortex solutions for a given
solution to the linear eigenvalue problem.

## Determining Solutions

We can calculate the streamfunction and vorticity fields in each layer given a
set of coefficients. These are determined by first calculating:

F = -U/ℓ sin θ sum_{j=0}^{M-1} a_j R_j(r/ℓ),

where (U, ℓ) denotes the vortex speed and radius, a_j is the j^th row of the
coefficient matrix `a`, M is the number of coefficients calculated and (r, θ)
are plane polar coordinates about the vortex center. The function R_j denotes
the Zernike radial function of order j.

The streamfunction and vorticity can be calculated from F using:

Δ_N(β) ψ = F,

and

q = Δ_N(0) ψ,

where

Δ_N(β) = [∇²-R[1]⁻²-β[1]/U,      R[1]⁻²     ,   0   , ...             ...        0;
	       R[2]⁻²     , ∇²-R[1]⁻²-β[1]/U, R[2]⁻², ...             ...        0;
	        ...       ,        ...      ,  ...  , ...   ... ,     ...         ;
	         0        ,         0       ,   0   , ... R[N]⁻², ∇²-R[N]⁻²-β[N]/U]

## Fourier Transforms

We calculate ψ and q in Fourier space where ∇² = -(k² + l²) using the FFTW
package. Since ψ and q are real, we use the `rfft` and `irfft` functions and
create wavenumber grids which are consistent with this choice.

This package contains a function for creating the grid (`CreateGrid`) however
these functions are designed to be consistent with the `TwoDGrid` from
FourierFlows.jl and GeophysicalFlows.jl so will also work with this grid.

"""

"""
    Calc_ψq(grid, a; U, ℓ, R, β, x₀=[0, 0], α=0)

Calculate ``ψ`` and ``q`` in a layered QG model using coefficients and vortex parameters

# Arguments:
 - `grid`: grid structure containing x, y, and Krsq
 - `a`: M x N array of coefficients, Array
 - (`U`, `ℓ`): vortex speed and radius, Numbers (default: `(1, 1)`)
 - (`R`, `β`): Rossby radii and (y) PV gradients in each layer, Numbers or Vectors (default: `(Inf, 0)`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)
 - `α`: initial angle of vortex, Number (default: `0`)
"""
function Calc_ψq(
    grid,
    a::Array;
    U::Number = 1,
    ℓ::Number = 1,
    R::Union{Number,Vector} = Inf,
    β::Union{Number,Vector} = 0,
    x₀::Vector = [0, 0],
    α::Number = 0,
)

    M, N = size(a)

    Nx, Ny = length(grid.x), length(grid.y)

    # Create Cartesian and polar grids

    x, y = CartesianGrid(grid)
    r, θ = PolarGrid(x, y, x₀)

    # Define temporary variable for RHS terms

    F = zeros(Nx, Ny, N)

    # Set RHS terms as sum of coefficients multiplied by Zernike polynomials

    for j = 1:M

        for i = 1:N

            F[:, :, i] += a[j, i] * ZernikeR(j - 1, r / ℓ)

        end

    end

    # Fourier transform F, set (0, 0) mode to 0 to avoid NaN errors later

    F = -U / ℓ * F .* sin.(θ .- α)
    Fh = rfft(F, [1, 2])
    Fh[1, 1, :] .= 0

    # Define PV inversion operators for ψ and q

    ΔN_i = stack(inv, eachslice(ΔNCalc(grid.Krsq, R, β, U), dims = (3, 4)))
    ΔN_0 = ΔNCalc(grid.Krsq, R, 0)

    # Define temporary variables for ψ and q in Fourier space

    if grid.Krsq isa CuArray

        Fh = CuArray(Fh)
        ψh = CuArray{ComplexF64}(zeros(Int(Nx / 2 + 1), Ny, N))
        qh = CuArray{ComplexF64}(zeros(Int(Nx / 2 + 1), Ny, N))

    else

        ψh = Array{ComplexF64}(zeros(Int(Nx / 2 + 1), Ny, N))
        qh = Array{ComplexF64}(zeros(Int(Nx / 2 + 1), Ny, N))

    end

    # Calculate ψ in each layer in Fourier space

    for n = 1:N

        for j = 1:N

            ψh[:, :, n] .+= ΔN_i[n, j, :, :] .* Fh[:, :, j]

        end

    end

    # Calculate q in each layer in Fourier space

    for n = 1:N

        for j = 1:N

            qh[:, :, n] .+= ΔN_0[n, j, :, :] .* ψh[:, :, j]

        end

    end

    # Transform back to real space

    ψ = irfft(ψh, Nx, [1, 2])
    q = irfft(qh, Nx, [1, 2])

    return ψ, q

end

"""
    ΔNCalc(K², R, β, U=1)

Defines the ``Δ_N(β)`` matrix used to invert for ``ψ`` and ``q``

# Arguments:
 - `K²`: value of ``k²+l²`` in Fourier space, Array
 - (`R`, `β`): Rossby radii and (y) PV gradients in each layer, Numbers or Vectors
 - `U`: vortex speed, Number (default: `1`)
"""
function ΔNCalc(
    K²::Union{CuArray,Array},
    R::Union{Number,Vector},
    β::Union{Number,Vector},
    U::Number = 1,
)

    N = length(R)
    Nk, Nl = size(K²)
    ϵ = max(1, maximum(R .^ -2)) * 1e-15

    # Define βU⁻¹ depending on input type

    if length(β) < N

        βU⁻¹ = β[1] / U * ones(N)

    else

        βU⁻¹ = β / U

    end

    if N == 1

        # Calculate ΔN in 1-layer case

        ΔN = -reshape(K² .+ (1 / R^2 + βU⁻¹ + ϵ), 1, 1, Nk, Nl)

    else

        # Calculate ΔN in N-layer case (N > 1)

        diagonal_elements = -[R[1]^-2; 2 * R[2:end-1] .^ -2; R[end]^-2] - βU⁻¹
        above_diagonal_elements = R[1:end-1] .^ -2
        below_diagonal_elements = R[2:end] .^ -2

        if K² isa CuArray

            ΔN = CuArray(zeros(N, N, Nk, Nl))

            K2 = reshape(K² .+ ϵ, 1, 1, Nk, Nl)

            [ΔN[i, i, :, :] = -K2 for i in range(1, N)]

            ΔN .+= CuArray(
                diagm(
                    0 => diagonal_elements,
                    1 => above_diagonal_elements,
                    -1 => below_diagonal_elements,
                ),
            )

        else

            ΔN =
                diagm(
                    0 => diagonal_elements,
                    1 => above_diagonal_elements,
                    -1 => below_diagonal_elements,
                ) .- I(N) .* reshape(K² .+ ϵ, 1, 1, Nk, Nl)

        end

    end

    return ΔN

end

"""
# Low-Level Functions

This file contains functions for calculating SQG vortex solutions for a given
solution to the linear eigenvalue problem.

## Determining Solutions

We calculate the streamfunction and vorticity fields in each layer given a set
of coefficients. These are determined by first calculating:

F = -U/ℓ sin θ sum_{j=0}^{M-1} a_j R_j(r/ℓ),

where (U, ℓ) denotes the vortex speed and radius, a_j is the j^th row of the
coefficient matrix `a`, M is the number of coefficients calculated and (r, θ)
are plane polar coordinates about the vortex center. The function R_j denotes
the Zernike radial function of order j.

The streamfunction and buoyancy can be calculated from F using:

[∂/∂z + 1/R'] ψ = -ℓF,

and

b = ∂ψ/∂z,

where ∂/∂z = √(-∇² + β/U) tanh[√(-∇² + β/U) R] is a Dirichlet-Neumann operator
which relates the surface derivative of ψ to it's surface value.

## Fourier Transforms

We calculate ψ and b in Fourier space where ∇² = -(k² + l²) using the FFTW
package. Since ψ and b are real, we use the `rfft` and `irfft` functions and
create wavenumber grids which are consistent with this choice.

This package contains a function for creating the grid (`CreateGrid`) however
these functions are designed to be consistent with the `TwoDGrid` from
FourierFlows.jl and GeophysicalFlows.jl so will also work with this grid.

"""

"""
    Calc_ψb(grid, a; U, ℓ, R, β, x₀=[0, 0], α=0)

Calculate SQG fields ``ψ`` and ``b`` using coefficients and vortex parameters

# Arguments:
 - `grid`: grid structure containing `x`, `y`, and `Krsq`
 - `a`: M x 1 array of coefficients, Array

# Keyword arguments:
 - (`U`, `ℓ`): vortex speed and radius, Numbers (default: `1`)
 - `R`: vector of ``[R, R']``, Vector (default: `[Inf, Inf]`)
 - `β`: beta-plane (y) PV gradient, Number (default: `1`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)
 - `α`: initial angle of vortex, Number (default: `0`)

Note: Here R is the baroclinic Rossby radius, R = NH/f, and R' = R₀²/R where R₀ is
the barotropic Rossby radius, R₀ = √(gH)/f. For infinite depth, R' = g/(fN).
"""
function Calc_ψb(
    grid,
    a::Array;
    U::Number = 1,
    ℓ::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
    x₀::Vector = [0, 0],
    α::Number = 0,
)

    M, _ = size(a)
    Nx, Ny = length(grid.x), length(grid.y)
    ϵ = 1e-15

    # Create Cartesian and polar grids

    x, y = CartesianGrid(grid)
    r, θ = PolarGrid(x, y, x₀)

    # Define temporary variable for RHS terms

    F = zeros(Nx, Ny)

    # Set RHS terms as sum of coefficients multiplied by Zernike polynomials

    for j = 1:M

        F[:, :] += a[j, 1] * ZernikeR(j - 1, r / ℓ)

    end

    # Fourier transform F, set (0, 0) mode to 0 to avoid NaN errors later

    F = U * F .* sin.(θ .- α)
    Fh = rfft(F)
    Fh[1, 1] = 0

    # Define Dirichlet-Neumann operator linking ψ and b

    ∂z = sqrt.(grid.Krsq .+ β / U) .* tanh.((ϵ .+ sqrt.(grid.Krsq .+ β / U)) .* R[1])

    # Move Fh to GPU if `cuda=true` in grid

    if grid.Krsq isa CuArray

        Fh = CuArray(Fh)

    end

    # Calculate ψ and b in Fourier using Dirichlet-Neumann operator

    ψh = Fh ./ (∂z .+ (1 / R[2] + ϵ))
    bh = ∂z .* ψh

    # Transform back to real space

    ψ = irfft(ψh, Nx)
    b = irfft(bh, Nx)

    return ψ, b

end

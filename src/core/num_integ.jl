"""
This file contains numerical integration functions:

* `JJ_int`: calculates the integrals in the linear systems for vortex solutions of
various QG problems. Functions specific to a particular model are included in files
related to that model.

* `AreaInteg2`: calculates the area integral of a function squared. Used for energy
and enstrophy diagnostics for various QG models.

"""

"""
    JJ_int(F, j, k, tol=1e-6)

Evaluates the integral ``I = ∫_0^∞ F(ξ) J_{2j+2}(ξ) J_{2k+2}(ξ) \\mathrm{d}ξ``

# Arguments:
 - `F`: function to integrate, typically `A_func` or `B_func`, Function
 - `j`: first Bessel function index, Integer
 - `k`: second Bessel function index, Integer
 - `tol`: error tolerance for QuadGK, Number (default: `1e-6`)

Note: This integral is performed by deforming the contour of integration into the complex plane
where the Bessel function decays exponentially in the imaginary direction.
"""
function JJ_int(F::Function, j::Int, k::Int, tol::Number = 1e-6)

    # Define parameters for numerical integration

    d, D = 1e3, 100# splitting parameter and domain limit for exp term
    atol = 1e-4 * tol# absolute error tolerance for quadgk

    # Define Bessel functions

    J(n, x) = besselj(n, x)
    Y(n, x) = bessely(n, x)
    H₁(n, x) = hankelh1(n, x)
    H₂(n, x) = hankelh2(n, x)

    # Define products of Bessel functions

    HH₁₂(m, n, x) = @. 2 * (J(m, x) * J(n, x) + Y(m, x) * Y(n, x))
    ϕ₁₁(m, n, x) = @. H₁(m, x) * H₁(n, x) / exp(2 * im * x)
    ϕ₂₂(m, n, x) = @. H₂(m, x) * H₂(n, x) / exp(-2 * im * x)

    # Define integrands for each of the 4 contour segments

    F₁(x) = F(x) .* J(2 * j + 2, x) * J(2 * k + 2, x)
    F₂(x) = F(x) .* HH₁₂(2 * j + 2, 2 * k + 2, x) / 4
    F₃(x) = F(d + im * x) .* (ϕ₁₁(2 * j + 2, 2 * k + 2, d + im * x) * exp(-2 * x))
    F₄(x) = F(d - im * x) .* (ϕ₂₂(2 * j + 2, 2 * k + 2, d - im * x) * exp(-2 * x))

    # Integrate along each contour segment

    I₁ = quadgk(F₁, 0, d, rtol = tol, atol = atol)
    I₂ = quadgk(F₂, d, Inf, rtol = tol, atol = atol)
    I₃ = quadgk(F₃, 0, D, rtol = tol, atol = atol)
    I₄ = quadgk(F₄, 0, D, rtol = tol, atol = atol)

    # Define integral I and total error ϵ

    I = I₁[1] + I₂[1] - imag(exp(2 * im * d) * I₃[1] - exp(-2 * im * d) * I₄[1]) / 4
    ϵ = I₁[2] + I₂[2] + I₃[2] + I₄[2]

    return I, ϵ

end

"""
    AreaInteg2(f, grid)

Calculates the integral ``I = ∫_A f^2 \\mathrm{d}A`` where ``A``
is the 2D domain described by `grid`.

# Arguments:
 - `f`: input Array in real or Fourier space
 - `grid`: grid structure

Note: f can be entered in real space or Fourier space, we use the rfft function
to calculate the Fourier transform so array sizes can distinguish the two.
"""
function AreaInteg2(f::Union{CuArray,Array}, grid::GridStruct)

    # Get grid parameters

    Nkr = length(grid.kr)
    Nx, Ny = length(grid.x), length(grid.y)
    Δx, Δy = grid.x[2] - grid.x[1], grid.y[2] - grid.y[1]

    # If input array is in real space, apply Fourier transform

    if size(f)[1:2] == (Nx, Ny)

        fh = rfft(f, [1, 2])

    else

        fh = f

    end

    # Add up components in Fourier space, non-edge elements are counted twice for an rfft

    I =
        sum(abs2, fh[1, :]) +             # kr = 0 (edge)
        sum(abs2, fh[Nkr, :]) +         # kr = Nx/2 (edge)
        2 * sum(abs2, fh[2:Nkr-1, :])  # sum twice for non-edge modes (as using rfft)

    I = I * (Δx * Δy) / (Nx * Ny) # normalization factor for fft

    return I

end

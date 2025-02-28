"""
This file contains the numerical integration functions required to calculate the terms A_{k,j} and
B_{k,j} for building the matrices A and B for the layered QG model. These terms are given by:

A_{k,j} = int_0^inf K(ξ) [K(ξ) + D(μ)]⁻¹ ξ⁻¹ J_{2j+2}(ξ) J_{2k+2}(ξ) dξ,

B_{k,j} = int_0^inf [K(ξ) + D(μ)]⁻¹ ξ⁻¹ J_{2j+2}(ξ) J_{2k+2}(ξ) dξ,

where J_n(ξ) denotes the Bessel function of order n, K(ξ) is given by

K(ξ) = [ξ²+λ[1]²,  -λ[1]²  ,   0   , ...            ... 0;
	 -λ[2]² , ξ²+2λ[2]², -λ[2]², ...            ... 0;
	  ...   ,    ...   ,  ...  , ...   ... ,   ...   ;
	    0   ,     0    ,   0   , ... -λ[N]², ξ²+λ[N]²]

and D(μ) = diag(μ[1], μ[2], ... μ[N]).

For the SQG model, the expressions are more straightforward and defined in `lin_sys.jl` when
required for building A and B. The function `JJ_int` below is still used in the SQG case.
"""


"""
    A_func(ξ, λ, μ)

Evaluates the matrix function ``A(ξ, λ, μ) = K(ξ) [K(ξ) + D(μ)]⁻¹ ξ⁻¹``

# Arguments:
 - `ξ`: point in ``[0, ∞)``, Number
 - `λ`: ratio of vortex radius to Rossby radius in each layer, Number or Vector
 - `μ`: nondimensional (y) vorticity gradient in each layer, Number or Vector
"""
function A_func(ξ::Number, λ::Union{Vector,Number}, μ::Union{Vector,Number})

    N = length(μ)

    if N == 1

        # Calculate A in 1-layer case

        K = @. ξ^2 + λ^2

        A = @. K / (K + μ) / ξ

    else

        # Calculate A in N-layer case (N > 1)

        diagonal_elements = [λ[1]^2; 2 * λ[2:end-1] .^ 2; λ[end]^2]
        above_diagonal_elements = -λ[1:end-1] .^ 2
        below_diagonal_elements = -λ[2:end] .^ 2

        K =
            I(N) * ξ^2 + diagm(
                0 => diagonal_elements,
                1 => above_diagonal_elements,
                -1 => below_diagonal_elements,
            )

        A = (K / (K .+ diagm(μ))) / ξ

    end

    return A

end

"""
    B_func(ξ, λ, μ)

Evaluates the matrix function ``B(ξ, λ, μ) = [K(ξ) + D(μ)]⁻¹ ξ⁻¹``

# Arguments:
 - `ξ`: point in ``[0, ∞)``, Number
 - `λ`: ratio of vortex radius to Rossby radius in each layer, Number or Vector
 - `μ`: nondimensional (y) vorticity gradient in each layer, Number or Vector
"""
function B_func(ξ::Number, λ::Union{Vector,Number}, μ::Union{Vector,Number})

    N = length(μ)

    if N == 1

        # Calculate B in 1-layer case

        K = @. ξ^2 + λ^2

        B = @. 1 / (K + μ) / ξ

    else

        # Calculate B in N-layer case (N > 1)

        diagonal_elements = [λ[1]^2; 2 * λ[2:end-1] .^ 2; λ[end]^2]
        above_diagonal_elements = -λ[1:end-1] .^ 2
        below_diagonal_elements = -λ[2:end] .^ 2

        K =
            I(N) * ξ^2 + diagm(
                0 => diagonal_elements,
                1 => above_diagonal_elements,
                -1 => below_diagonal_elements,
            )

        B = inv(K .+ diagm(μ)) / ξ

    end

    return B

end

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

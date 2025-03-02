"""
# Low-Level Functionality for SQG Model

## The Linear System:

For the SQG problem, A_{j,k} = δ_{j,k} / (4*j), B is determined using `JJ_int` as a double Bessel
integral of F(ξ) = [D(ξ) ξ]⁻¹ where:

D(ξ) = sqrt(ξ^2 + μ) * tanh(sqrt(ξ^2 + μ) / λ[1]) + λ[2], 	for λ[1] > 0
       sqrt(ξ^2 + μ) + λ[2],					for λ[1] = 0

c_i = 1/4 δ_{i, 1},

c₀ = 0,

and

d_i = (-1)^i.

Note that for λ = [0, 0] and μ = 0, B can be calculated analytically as:

B_{j+1, k+1} = 4*(-1)^(j-k+1)/((2j-2k-1)*(2j-2k+1)*(2j+2k+3)*(2j+2k+5))/π.

## Solving the System:

The linear system can be solved using `SolveInhomEVP` to obtain the eigenvalues/
eigenvectors.

## Calculating the Solution:

The solution may be calculated using functions from `sqg_low.jl`.

"""

"""
    BuildLinSysSQG(M, λ, μ; tol=1e-6)

Builds the terms in the inhomogeneous eigenvalue problem; ``A``, ``B``, ``c``, ``d`` for
the SQG problem

# Arguments:
 - `M`: number of coefficient to solve for, Integer
 - `λ`: ratio of vortex radius to Rossby radius in each layer, Number or Vector
 - `μ`: nondimensional (y) vorticity gradient in each layer, Number or Vector
 - `tol`: error tolerance for QuadGK via `JJ_int`, Number (default: `1e-6`)
"""
function BuildLinSysSQG(
    M::Int,
    λ::Union{Vector,Number},
    μ::Union{Vector,Number};
    tol::Number = 1e-6,
)

    # Define A, B, c, d for the SQG system, B is temporary here

    A = diagm(1 ./ (1:M) / 4)
    B = zeros(M, M)
    c = hcat(zeros(M, 1), vcat(1 / 4, zeros(M - 1, 1)))
    d = reshape((-1) .^ (0:M-1), M, 1)

    # Set values of B depending on cases

    if (μ == 0) & (λ == [0, 0])

        # When analytic solution exists, use that

        B₀(j, k) =
            4 * (-1)^(j - k + 1) /
            ((2j - 2k - 1) * (2j - 2k + 1) * (2j + 2k + 3) * (2j + 2k + 5)) / π

        for j = 0:M-1

            for k = 0:M-1

                B[j+1, k+1] = B₀(j, k)

            end

        end

    else

        # When analytic solution does not exists, use numerical integration

        D_func(ξ) = @. sqrt(ξ^2 + μ) * tanh(sqrt(ξ^2 + μ) / λ[1]) + λ[2]# SQG kernel function

        for j = 0:M-1

            for k = 0:M-1

                B[j+1, k+1] = JJ_int(x -> 1 ./ (D_func(x) .* x), j, k, tol)[1]
            end

        end

    end

    return A, B, c, d

end

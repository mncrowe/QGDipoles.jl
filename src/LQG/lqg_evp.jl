"""
# Low-Level Functionality for LQG Model

## The Linear System:

For the LQG problem, A and B₀ are given by the terms A_{j,k}, B_{j,k} and B[n] contains only 
the rows of B₀ corresponding to the coefficients in the n^th layer.  The vectors c₀, c[n] and
d[n] are given by:

(c[n])_i = 1/4 δ_{i, n},

c₀ = sum_{n=1}^N (μ[n] c[n]),

and

(d[n])_i = (-1)^(i/n) * δ_{mod(i, N), n}.

The terms A_{k,j} and B_{k,j} are given by:

A_{k,j} = int_0^inf K(ξ) [K(ξ) + D(μ)]⁻¹ ξ⁻¹ J_{2j+2}(ξ) J_{2k+2}(ξ) dξ,
B_{k,j} = int_0^inf [K(ξ) + D(μ)]⁻¹ ξ⁻¹ J_{2j+2}(ξ) J_{2k+2}(ξ) dξ,

where J_n(ξ) denotes the Bessel function of order n, K(ξ) is given by

K(ξ) = [ξ²+λ[1]²,  -λ[1]²  ,   0   , ...            ... 0;
	 -λ[2]² , ξ²+2λ[2]², -λ[2]², ...            ... 0;
	  ...   ,    ...   ,  ...  , ...   ... ,   ...   ;
	    0   ,     0    ,   0   , ... -λ[N]², ξ²+λ[N]²]

and D(μ) = diag(μ[1], μ[2], ... μ[N]).

## Active and Passive Layers:

Layers without a vortex core can be included in the model. They can be removed from the 
linear system using `ApplyPassiveLayers` and then included back into the eigenvalue/
eigenvector solution using `IncludePassiveLayers`.

## Solving the System:

The linear system can be solved using `SolveInhomEVP` to obtain the eigenvalues/
eigenvectors.

## Calculating the Solution:

The solution may be calculated using functions from `lqg_low.jl`.

"""

"""
    BuildLinSysLQG(M, λ, μ; tol=1e-6)

Builds the terms in the inhomogeneous eigenvalue problem; ``A``, ``B``, ``c``, ``d`` for
the LQG problem

# Arguments:
 - `M`: number of coefficient to solve for, Integer
 - `λ`: ratio of vortex radius to Rossby radius in each layer, Number or Vector
 - `μ`: nondimensional (y) vorticity gradient in each layer, Number or Vector

# Keyword arguments:
 - `tol`: error tolerance for QuadGK via `JJ_int`, Number (default: `1e-6`)
"""
function BuildLinSysLQG(
    M::Int,
    λ::Union{Vector,Number},
    μ::Union{Vector,Number};
    tol::Number = 1e-6,
)

    # Get number of layers from input size

    N = length(μ)

    # Set temporary values for A, B, c, d for LQG case

    A = zeros(N * M, N * M)
    B₀ = zeros(N * M, N * M)
    B = zeros(N * M, N * M, N)
    c = zeros(N * M, N + 1)
    d = zeros(N * M, N)
    c₀ = vcat(ones(N), zeros((M - 1) * N))

    # Set A and temporary B values using numerical integration

    for j = 0:M-1

        for k = 0:M-1

            A[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> A_func(x, λ, μ), j, k, tol)[1]
            B₀[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> B_func(x, λ, μ), j, k, tol)[1]

        end

    end

    # Set B, c, d values using existing quantities

    for n = 1:N

        K = kron(I(M), diagm((1:N) .== n))

        B[:, :, n] = K * B₀
        c[:, 1] = c[:, 1] + μ[n] * (K * c₀) / 4
        c[:, n+1] = (K * c₀) / 4
        d[:, n] = kron((-1) .^ (0:M-1), (1:N) .== n)

    end

    return A, B, c, d

end

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
    ApplyPassiveLayers(A, B, c, d, ActiveLayers)

Removes rows and columns corresponding to passive layers from the system

# Arguments:
 - `A`, `B`, `c`, `d`: inhomogeneous eigenvalue problem terms, Arrays
 - `ActiveLayers`: vector of 1s or 0s where 1 denotes an active layer, Number or Vector
"""
function ApplyPassiveLayers(
    A::Array,
    B::Array,
    c::Array,
    d::Array,
    ActiveLayers::Union{Number,Vector},
)

    # Ensure ActiveLayers is a Vector

    if ActiveLayers isa Number

        ActiveLayers = [ActiveLayers]

    end

    # Calculate number of coefficients from size of input

    M = Int(size(d)[1] / size(d)[2])# problem size

    # Define arrays of true/false flags for active and passive layers

    i₁ = BitArray{1}(kron(ones(M), ActiveLayers))# grid index of active layers
    i₃ = BitArray{1}(ActiveLayers)# index of active layers
    i₄ = BitArray{1}(vcat(1, ActiveLayers))# extended index of active layers

    # Define new system by removing rows/columns corresponding to passive layers

    A = A[i₁, i₁]
    B = B[i₁, i₁, i₃]
    c = c[i₁, i₄]
    d = d[i₁, i₃]

    return A, B, c, d

end

"""
    IncludePassiveLayers(K, a, ActiveLayers)

Includes columns corresponding to passive layers in the eigenvalue and coefficient arrays

# Arguments:
 - `K`, `a`: eigenvalue and coefficient arrays describing system solution, Arrays
 - `ActiveLayers`: vector of 1s or 0s where 1 denotes an active layer, Number or Vector
"""
function IncludePassiveLayers(K::Array, a::Array, ActiveLayers::Union{Number,Vector})

    # Ensure ActiveLayers is a Vector

    if ActiveLayers isa Number

        ActiveLayers = [ActiveLayers]

    end

    # Get numbers of coefficients and layers using input size

    M = size(a)[1]
    N = length(ActiveLayers)

    # Define variables for (K, a) corresponding to full system size

    K₁, a₁ = zeros(1, N), zeros(M, N)

    # Set (K, a) values in active layers, passive layers retain values of 0

    i = BitArray{1}(ActiveLayers)

    K₁[:, i] .= K
    a₁[:, i] .= a

    return K₁, a₁

end

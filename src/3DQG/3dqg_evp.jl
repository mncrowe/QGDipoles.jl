"""
This file contains functions for building the eigenvalue problem for the 3D QG model

...

"""




"""
    BuildLinSys3DQG(M, Nz, λ, μ, λ₀, γ; tol=1e-6)

Builds the terms in the inhomogeneous eigenvalue problem; ``A``, ``B``, ``c``, ``d`` for
the LQG problem

# Arguments:
 - `M`: number of coefficient to solve for, Integer
 - `N`: number of z points, Integer
 - `λ`: ratio of vortex radius to baroclinic Rossby radius at each height ζ, Function
 - `μ`: nondimensional (y) vorticity gradient in each layer, Number
 - `λ₀`: ratio of baroclinic Rossby radius to barotropic Rossby radius at surface, Number
 - `γ`: nondimensional (y) slope on bottom surface, Number

# Keyword arguments:
 - `tol`: error tolerance for `QuadGK` via `JJ_int`, Number (default: `1e-6`)
"""
function BuildLinSys3DQG(
    M::Int,
    N::Int,
    λ::Function,
    μ::Number,
    λ₀::Number,
    γ::Number;
    tol::Number = 1e-6,
)

    # Create ζ grid

    ζ, D_ζ = GridChebyshev(N, [-1, 0])

    # Create matrices for evalating A and B

    μ_vec = [-γ; ones(N-2); 0]
    Λ = diagm(λ.(ζ) .^ 2)

    M₀ = - D_ζ * Λ * D_ζ
    M₀[1, :] = D_ζ[1, :]
    M₀[N, :] = (D_ζ + λ₀^2 * I(N))[N, :]

    I₀ = diagm([0; ones(N-2); 0])
    C = diagm(μ_vec)

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
        c[:, 1] = c[:, 1] + μ_vec[n] * (K * c₀) / 4
        c[:, n+1] = (K * c₀) / 4
        d[:, n] = kron((-1) .^ (0:M-1), (1:N) .== n)

    end

    return A, B, c, d

end



function A_func(ξ::Number, M₀::Matrix, I₀::Matrix, C::Matrix)

    M_ξ = ξ^2 * I₀ + M₀

    A = M_ξ / (M_ξ + C) / ξ

    return A

end




function B_func(ξ::Number, M₀::Matrix, I₀::Matrix, C::Matrix)

    M_ξ = ξ^2 * I₀ + M₀

    B = inv(M_ξ + C) / ξ

    return A

end




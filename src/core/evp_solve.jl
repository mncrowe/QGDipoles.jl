"""
This file contains the numerical integration functions required to build the inhomogeneous eigenvalue
problem. This problem is of size MN x MN where M is the number of coefficients in each layer and N is
the number of layers. The problem is:

    [A - sum_{n=1}^N (Kᵐ[n] B[n])] a = c₀ + sum_{n=1}^N (Kᵐ[n] c[n]),  d[j]ᵀa = 0, j ∈ {1, .., N}

where eigenvalues are denoted by Kᵐ[n] and the eigenvector by `a`. The exponent m is determined by the
problem type; the layered QG model has m = 2 while the SQG problem has m = 1.

This system is solved using two approaches. For N = 1 (inc. SQG problem), the system may be
converted into a generalised eigenvalue problem of size 2M x 2M and solved directly. For N > 1,
we instead use a nonlinear root finding approach with the NLSolve package. The method works by
projecting `a` onto the subspace perpendicular to the d[n] vectors. A vector x is defined as
x = [Kᵐ; a'] where a' denotes the projection of a. Since a' has N less degrees of freedom than
a, and K² is of length N, the vector x is of length M*N. Defining

F(x) = [A - sum_{n=1}^N (Kᵐ[n] B[n])] a - c₀ - sum_{n=1}^N (Kᵐ[n] c[n]),

allows us to solve the inhomogeneous problem by finding roots of F(x) = 0 using some initial guess
x₀ = [Kᵐ₀; a₀']. Changing the initial guess may be required to identify the required solutions.

"""

"""
    SolveInhomEVP(A, B, c, d; K₀=nothing, a₀=nothing, tol=1e-6, method=0, m=2, warn=true)

Solves the inhomogeneous eigenvalue problem using nonlinear root finding

# Arguments:
 - `A`, `B`, `c`, `d`: inhomogeneous eigenvalue problem terms, Arrays
 - `K₀`, `a₀`: initial guesses for ``K`` and ``a``, Arrays or nothings (default: `nothing`)
 - `tol`: error tolerance for `nlsolve`, Number (default: `1e-6`)
 - `method`: `0` - eigensolve for ``N = 1`` and `nlsolve` for ``N > 1``, `1` - `nlsolve` (default: `0`)
 - `m`: exponent of ``K`` in eignevalue problem (default: `2`)
 - `warn`: if `true` displays warning if solution includes unextracted passive layers (default: `true`)

"""
function SolveInhomEVP(
    A::Array,
    B::Array,
    c::Array,
    d::Array;
    K₀::Union{Number,Array,Nothing} = nothing,
    a₀::Union{Array,Nothing} = nothing,
    tol::Number = 1e-6,
    method::Int = 0,
    m::Int = 2,
    warn::Bool = true,
)

    # Ensure K₀ is a Vector

    if K₀ isa Number

        K₀ = [K₀]

    end

    # Determine number of layers and number of coefficients from input

    N = size(d)[2]
    M = Int(size(d)[1] / N)

    # Use root finding method for multi-layer systems

    if N > 1

        method = 1

    end

    if method == 0

        # Use analytical solution if method = 0

        # Set K₀ value if none given

        if K₀ isa Nothing

            K₀ = [4]

        end

        # Reformat inputs as arrays with correct shape

        B = reshape(B, M, M)
        O = zeros(M, M)
        dᵀ = permutedims(d)
        c₀ = c[:, 1]
        c₁ = c[:, 2]

        # Build intermediate matrices for quadratic eigenvalue problems

        D₀ = (dᵀ * (A \ c₀)) .* A
        D₁ = (dᵀ * (A \ c₁)) .* A - (dᵀ * (A \ c₀)) .* B + (c₀ * dᵀ) * (A \ B)
        D₂ = -(dᵀ * (A \ c₁)) .* B + (c₁ * dᵀ) * (A \ B)

        # Convert quadratic system to canonical form

        D₃ = [D₀ O; O I(M)]
        D₄ = [-D₁ -D₂; I(M) O]

        # Solve for eigenvalues of quadratic eigenvalue problem

        λ = eigvals(D₃, D₄)

        # Identify eigenvalue closest to K₀^m

        v = abs.((λ .- K₀ .^ m) .^ 2)
        i = argmin(v[.!isnan.(v)])

        # Define K as requested eigenvalue and invert system for a

        K = reshape([λ[i]], 1, 1) .^ (1 / m)
        a = reshape((A - K .^ m .* B) \ (c₀ + K .^ m .* c₁), M, 1)

    end

    if method == 1

        # Use root finding method if method = 1

        # Calculate basis vectors spanning the space orthogonal to the d vectors

        e, V, iₑ = OrthogSpace(d)

        # Define a₀ if none given and reshape if given

        if a₀ isa Nothing

            a₀ = vcat(-10 * ones(N, 1), zeros(N * (M - 1), 1))

        else

            a₀ = reshape(permutedims(a₀), N * M, 1)

        end

        # Define K₀ if none given and reshape if given

        if K₀ isa Nothing

            K₀ = 5 * ones(N, 1)

        else

            K₀ = reshape(K₀, N, 1)

        end

        # Calculate x₀ by projecting a₀ onto the new space

        x₀ = V \ a₀

        # Define x₀ as vector of K values and a₀ values in the new space

        x₀ = vcat(K₀ .^ m, x₀[iₑ])

        # Define function and derivative for root finding

        fj! = (F, J, x) -> InhomEVP_F!(F, J, x, A, B, c, e)

        # Solve for F(x) = 0

        x = nlsolve(only_fj!(fj!), x₀, ftol = tol).zero

        # Extract K values from first few entries of x

        K = (complex(reshape(x[1:N], 1, N))) .^ (1 / m)

        # Calculate a values by projecting x back to original space

        a = permutedims(reshape(e * x[N+1:N*M], N, M))

    end

    # Raise warning if root finding has converged to (possibly) wrong solution unless suppressed

    if (imag(K) != zeros(1, N)) & warn

        @warn "Solution has complex K, generally corresponding passive layers."

    end

    # Remove any imaginary parts of K and ignore coefficients less than `tol` value

    K = real(K)
    a[abs.(a).<tol] .= 0

    return K, a

end

"""
    InhomEVP_F!(F, J, x, A, B, c, d, e)

Calculates the function ``F`` and it's derivatives, ``J``, at a given point ``x``

# Arguments:
 - `F`, `J`: values of ``F`` and it's derivatives, updated by function
 - `x`: evaluation point, Array
 - `A`, `B`, `c`: inhomogeneous eigenvalue problem terms, Arrays
 - `e`: basis spanning the space perpendicular to the ``d[n]``, Array
"""
function InhomEVP_F!(F, J, x::Array, A::Array, B::Array, c::Array, e::Array)

    # Get problem size from inputs

    N, j = size(e)

    # Calculate a by projecting x into original space

    a = e * x[N-j+1:N]

    # Build matrix and RHS vector for F(x)

    M = A
    v = c[:, 1]

    for i = 1:N-j

        M = M - x[i] * B[:, :, i]
        v = v + x[i] * c[:, i+1]

    end

    # Calculate and update Jacobian matrix of derivatives

    if !(J == nothing)

        for i = 1:N-j

            J[:, i] = -B[:, :, i] * a - c[:, i+1]

        end

        J[:, N-j+1:end] .= M * e

    end

    # Calculate and update function value

    if !(F == nothing)

        F[:] .= M * a - v

    end

end

"""
    OrthogSpace(v)

Extends the input to an orthonormal basis over ``R^n`` using the Gram-Schmidt method

# Arguments:
 - `v`: array with vectors as columns, Array
"""
function OrthogSpace(v)

    # Get problem size from size of imputs

    N = size(v)[1]

    if length(size(v)) > 1

        k = size(v)[2]

    else

        k = 1

    end

    # Set orthononality threshold for determining linear independence

    ϵ = 1e-6

    # Set temporary values of basis vectors, B, and new vector indices

    B = Matrix{Float64}(I, N, N)
    iₑ = 1:N

    # Perform the Gram-Schmidt algorithm to extend v to a full basis of R^N

    for i = 1:k

        j = 1

        while length(iₑ) > N - i

            if j > length(iₑ)

                @error "The v must be linerly independent."

            end

            if dot(v[:, i], B[:, iₑ[j]]) > ϵ

                B[:, iₑ[j]] = v[:, i]
                iₑ = setdiff(iₑ, iₑ[j])

            end

            j = j + 1

        end

    end

    # Orthogonalise full basis

    for j = 1:N

        for i = 1:j-1

            B[:, j] = B[:, j] - B[:, i] * dot(B[:, i], B[:, j]) / norm(B[:, i])^2

        end

    end

    # Normalise all vectors

    B = B ./ sqrt.(sum(abs2, B, dims = 1))

    # Define e as the new vectors added to the basis (B is the full basis)

    e = B[:, iₑ]

    return e, B, iₑ

end
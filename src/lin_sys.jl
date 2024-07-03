"""
This file contains the numerical integration functions required to build the inhomogeneous eigenvalue
problem:

    [A - sum_{n=1}^N (K²[n] B[n])] a = c₀ + sum_{n=1}^N (K²[n] c[n]),  d[j]ᵀa = 0, j ∈ {1, .., N}

A and B₀ are given by the terms A_{j,k}, B_{j,k} in `JJ_int.jl` and B[n] contains only the rows of B₀
corresponding to the coefficients in the n^th layer. The eigenvalues are denoted by K²[n] and the
eigenvector by `a`. The vectors c₀, c[n] and d[n] are given by:

(c[n])_i = 1/4 δ_{i, n},

c₀ = sum_{n=1}^N (μ[n] c[n]),

and

(d[n])_i = (-1)^(i/2) * δ_{mod(i, N), n}.

This system is solved with nonlinear root finding using the NLSolve package. The method works by
projecting `a` onto the subspace perpendicular to the d[n] vectors. A vector x is defined as
x = [K²; a'] where a' denotes the projection of a. Since a' has N less degrees of freedom than
a, and K² is of length N, the vector x is of length M*N where M is the number of coefficients in
each layer and N is the number of layers. Defining

F(x) = [A - sum_{n=1}^N (K²[n] B[n])] a - c₀ - sum_{n=1}^N (K²[n] c[n]),

allows us to solve the inhomogeneous problem by finding roots of F(x) = 0 using some initial guess
x₀ = [K²₀; a₀']. Changing the initial guess may be required to identify the required solutions.
"""


"""
Function: BuildLinSys

Builds the terms in the inhomogeneous eigenvalue problem; A, B, c, d

Arguments:
 - M: number of coefficient to solve for, Integer
 - λ: ratio of vortex radius to Rossby radius in each layer, Number or Vector
 - μ: nondimensional (y) vorticity gradient in each layer, Number or Vector
 - tol: error tolerance for QuadGK via `JJ_int`, Number (default: 1e-6)
"""

function BuildLinSys(M::Int, λ::Union{Vector,Number}, μ::Union{Vector,Number}; tol::Number=1e-6)

	N = length(μ)
	A, B₀ = zeros(N*M, N*M), zeros(N*M, N*M)

	for j in 0:M-1
		for k in 0:M-1

			A[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> A_func(x, λ, μ), j, k, tol)[1]
			B₀[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> B_func(x, λ, μ), j, k, tol)[1]

		end
	end

	B, c, d = zeros(N*M, N*M, N), zeros(N*M, N+1), zeros(N*M, N)
	c₀ = vcat(ones(N), zeros((M-1)*N))
	
	for n in 1:N
		
		K = kron(I(M), diagm((1:N).==n))
		B[:, :, n] = K * B₀
		c[:, 1] = c[:, 1] + μ[n] * (K * c₀) / 4
		c[:, n+1] = (K * c₀) / 4
		d[:, n] = kron((-1).^(0:M-1), (1:N).==n)
		
	end

	return A, B, c, d

end

"""
Function: ApplyPassiveLayers

Removes rows and columns corresponding to passive layers from the system

Arguments:
 - A, B, c, d: inhomogeneous eigenvalue problem terms, Arrays
 - ActiveLayers: vector of 1s or 0s where 1 denotes an active layer, Vector
"""

function ApplyPassiveLayers(A::Array, B::Array, c::Array, d::Array, ActiveLayers::Vector)

	M = Int(size(d)[1]/size(d)[2])			# problem size

	i₁ = BitArray{1}(kron(ones(M), ActiveLayers))	# grid index of active layers
	i₂ = BitArray{1}(1 .- ActiveLayers)		# index of passive layers
	i₃ = BitArray{1}(ActiveLayers)			# index of active layers
	i₄ = BitArray{1}(vcat(1, ActiveLayers))		# extended index of active layers
	
	A = A[i₁, i₁]
	B = B[i₁, i₁, i₃]
	c = c[i₁, i₄]
	d = d[i₁, i₃]
	
	return A, B, c, d
	
end

"""
Function: IncludePassiveLayers

Includes columns corresponding to passive layers in the eigenvalue and coefficient arrays

Arguments:
 - K, a: eigenvalue and coefficient arrays describing system solution, Arrays
 - ActiveLayers: vector of 1s or 0s where 1 denotes an active layer, Vector
"""

function IncludePassiveLayers(K::Array, a::Array, ActiveLayers::Vector)
	
	M, N = size(a)[1], length(ActiveLayers)

	K₁, a₁ = zeros(1, N), zeros(M, N)
	
	i = BitArray{1}(ActiveLayers)

	K₁[:, i] .= K
	a₁[:, i] .= a
	
	return K₁, a₁

end

"""
Function: SolveInhomEVP

Solves the inhomogeneous eigenvalue problem using nonlinear root finding

Arguments:
 - A, B, c, d: inhomogeneous eigenvalue problem terms, Arrays
 - K₀, a₀: initial guesses for K and a, Arrays or Nothings
 - tol: error tolerance for NLSolve, Number
"""

function SolveInhomEVP(A::Array, B::Array, c::Array, d::Array; K₀=Nothing, a₀=Nothing, tol::Number=1e-6)

	e, V, iₑ = OrthogSpace(d)

	N = size(d)[2]
	M = Int(size(d)[1]/N)
	
	if a₀ == Nothing
		a₀ = vcat(-10*ones(N, 1), zeros(N*(M-1), 1))
	else
		a₀ = reshape(permutedims(a₀), N*M, 1)
	end

	if K₀ == Nothing
		K₀ = 5*ones(N, 1)
	else
		K₀ = reshape(K₀, N, 1)
	end

	x₀ = V \ a₀
	x₀ = vcat(K₀.^2, x₀[iₑ])

	fj! = (F, J, x) -> InhomEVP_F!(F, J, x, A, B, c, e)

	x = nlsolve(only_fj!(fj!), x₀, ftol=tol).zero

	K = sqrt.(complex(reshape(x[1:N], 1, N)))
	a = permutedims(reshape(e * x[N+1:N*M], N, M))
	
	if imag(K) != zeros(1, N)
		@warn "Solution contains passive layers."
	end

	K = real(K)
	a[abs.(a) .< tol] .= 0

	return K, a

end

"""
Function: InhomEVP_F!

Calculates the function F and it's derivatives, J, at a given point x

Arguments:
 - F, J: values of F and it's derivatives, updated by function
 - x: evaluation point, Array
 - A, B, c: inhomogeneous eigenvalue problem terms, Arrays
 - e: basis spanning the space perpendicular to the d[n], Array
"""

function InhomEVP_F!(F, J, x::Array, A::Array, B::Array, c::Array, e::Array)

	N, j = size(e)

	a = e * x[N-j+1:N]
	M = A
	v = c[:, 1]
		
	for i in 1:N-j
		M = M - x[i] * B[:, :, i]
		v = v + x[i] * c[:, i + 1]
	end

	if !(J == nothing)

		for i in 1:N-j
			J[:, i] = -B[:, :, i] * a - c[:, i + 1]
		end

		J[:, N-j+1:end] .= M * e
		
	end
	
	if !(F == nothing)
		
		F[:] .= M * a - v

	end

end

"""
Function: OrthogSpace

Extends the input to an orthonormal basis over R^n using the Gram-Schmidt method

Arguments:
 - v: array with vectors as columns, Array
"""

function OrthogSpace(v)
	
	N = size(v)[1]
	if length(size(v)) > 1
		k = size(v)[2]
	else
		k = 1
	end

	ϵ = 1e-6

	B = Matrix{Float64}(I, N, N)
	iₑ = 1:N

	for i in 1:k
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

	for j in 1:N
		for i in 1:j-1
			B[:, j] = B[:, j] - B[:, i] * dot(B[:, i], B[:, j]) / norm(B[:, i])^2
		end
	end
	
	B = B ./ sqrt.(sum(abs2, B, dims=1))
	e = B[:, iₑ]
	
	return e, B, iₑ

end
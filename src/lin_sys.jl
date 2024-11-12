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

------------------------------------------------------------------------------------------------------

For the layered QG problem, A and B₀ are given by the terms A_{j,k}, B_{j,k} in `JJ_int.jl` and B[n]
contains only the rows of B₀ corresponding to the coefficients in the n^th layer.  The vectors c₀, 
c[n] and d[n] are given by:

(c[n])_i = 1/4 δ_{i, n},

c₀ = sum_{n=1}^N (μ[n] c[n]),

and

(d[n])_i = (-1)^(i/n) * δ_{mod(i, N), n}.

------------------------------------------------------------------------------------------------------

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
"""


"""
Function: `BuildLinSys(M, λ, μ; tol=1e-6, sqg=false)`

Builds the terms in the inhomogeneous eigenvalue problem; A, B, c, d

Arguments:
 - `M`: number of coefficient to solve for, Integer
 - `λ`: ratio of vortex radius to Rossby radius in each layer, Number or Vector
 - `μ`: nondimensional (y) vorticity gradient in each layer, Number or Vector
 - `tol`: error tolerance for QuadGK via `JJ_int`, Number (default: `1e-6`)
 - `sqg`: `false`; creates layered QG system, `true`; creates SQG system (default: `false`)
"""
function BuildLinSys(M::Int, λ::Union{Vector,Number}, μ::Union{Vector,Number}; tol::Number=1e-6, sqg::Bool=false)

	if sqg

		# Define A, B, c, d for the SQG system, B is temporary here

		A, B = diagm(1 ./(1:M)/4), zeros(M, M)

		c, d = hcat(zeros(M, 1), vcat(1/4, zeros(M-1, 1))), reshape((-1).^(0:M-1), M, 1)

		# Set values of B depending on cases
		
		if (μ == 0) & (λ == [0, 0])

			# When analytic solution exists, use that
			
			B₀(j, k) = 4*(-1)^(j-k+1)/((2j-2k-1)*(2j-2k+1)*(2j+2k+3)*(2j+2k+5))/π

			[[B[j+1, k+1] = B₀(j, k) for j = 0:M-1] for k = 0:M-1]
			
		else

			# When analytic solution does not exists, use numerical integration
			
			D_func(ξ) = @. sqrt(ξ^2 + μ) * tanh(sqrt(ξ^2 + μ) / λ[1]) + λ[2]

			[[B[j+1, k+1] = JJ_int(x -> 1 ./(D_func(x).*x), j, k, tol)[1] for j = 0:M-1] for k = 0:M-1]
			
		end

	else

		# Set values of A, B, c, d for LQG case

		N = length(μ)

		A, B₀ = zeros(N*M, N*M), zeros(N*M, N*M)

		# Set A and temporary B values using numerical integration

		[[A[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> A_func(x, λ, μ), j, k, tol)[1] for j = 0:M-1] for k = 0:M-1]

		[[B₀[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> B_func(x, λ, μ), j, k, tol)[1] for j = 0:M-1] for k = 0:M-1]

		# Define temporary values for B, c, d

		B, c, d = zeros(N*M, N*M, N), zeros(N*M, N+1), zeros(N*M, N)
		c₀ = vcat(ones(N), zeros((M-1)*N))

		# Set B, c, d values using existing quantities
	
		for n in 1:N
			
			K = kron(I(M), diagm((1:N).==n))

			B[:, :, n] = K * B₀

			c[:, 1] = c[:, 1] + μ[n] * (K * c₀) / 4
			c[:, n+1] = (K * c₀) / 4

			d[:, n] = kron((-1).^(0:M-1), (1:N).==n)
			
		end

	end

	return A, B, c, d

end

"""
Function: `ApplyPassiveLayers(A, B, c, d, ActiveLayers)`

Removes rows and columns corresponding to passive layers from the system

Arguments:
 - `A`, `B`, `c`, `d`: inhomogeneous eigenvalue problem terms, Arrays
 - `ActiveLayers`: vector of 1s or 0s where 1 denotes an active layer, Number or Vector
"""
function ApplyPassiveLayers(A::Array, B::Array, c::Array, d::Array, ActiveLayers::Union{Number,Vector})

	# Ensure ActiveLayers is a Vector

	if ActiveLayers isa Number
		
		ActiveLayers = [ActiveLayers]
		
	end

	# Calculate number of coefficients from size of input

	M = Int(size(d)[1]/size(d)[2])			# problem size

	# Define arrays of true/false flags for active and passive layers

	i₁ = BitArray{1}(kron(ones(M), ActiveLayers))	# grid index of active layers
	i₂ = BitArray{1}(1 .- ActiveLayers)		# index of passive layers
	i₃ = BitArray{1}(ActiveLayers)			# index of active layers
	i₄ = BitArray{1}(vcat(1, ActiveLayers))		# extended index of active layers
	
	# Define new system by removing rows/columns corresponding to passive layers

	A = A[i₁, i₁]
	B = B[i₁, i₁, i₃]
	c = c[i₁, i₄]
	d = d[i₁, i₃]
	
	return A, B, c, d
	
end

"""
Function: `IncludePassiveLayers(K, a, ActiveLayers)`

Includes columns corresponding to passive layers in the eigenvalue and coefficient arrays

Arguments:
 - `K`, `a`: eigenvalue and coefficient arrays describing system solution, Arrays
 - `ActiveLayers`: vector of 1s or 0s where 1 denotes an active layer, Number or Vector
"""
function IncludePassiveLayers(K::Array, a::Array, ActiveLayers::Union{Number,Vector})

	# Ensure ActiveLayers is a Vector

	if ActiveLayers isa Number
		
		ActiveLayers = [ActiveLayers]
		
	end

	# Get numbers of coefficients and layers using input size

	M, N = size(a)[1], length(ActiveLayers)

	# Define variables for (K, a) corresponding to full system size

	K₁, a₁ = zeros(1, N), zeros(M, N)
	
	# Set (K, a) values in active layers, passive layers retain values of 0

	i = BitArray{1}(ActiveLayers)

	K₁[:, i] .= K
	a₁[:, i] .= a
	
	return K₁, a₁

end

"""
Function: `SolveInhomEVP(A, B, c, d; K₀=Nothing, a₀=Nothing, tol=1e-6, method=0, m=2, sqg=false, warn=true)`

Solves the inhomogeneous eigenvalue problem using nonlinear root finding

Arguments:
 - `A`, `B`, `c`, `d`: inhomogeneous eigenvalue problem terms, Arrays
 - `K₀`, `a₀`: initial guesses for K and a, Arrays or Nothings (default: `Nothing`)
 - `tol`: error tolerance for `nlsolve`, Number (default: `1e-6`)
 - `method`: `0` - eigensolve for N = 1 and `nlsolve` for N > 1, `1` - `nlsolve` (default: `0`)
 - `m`: exponent of K in eignevalue problem (default: `2`)
 - `sqg`: `false`, uses `m` value specified; `true`, sets `m=1` (default: `false`)
 - `warn`: if `true` displays warning if solution includes unextracted passive layers (default: `true`)

Note: setting `sqg=true` overwrites the value of `m` and is equivalent to setting `m=1`.
The option to set both is included for consistency with `BuildLinSys` and more generality
with the value of `m`.
"""
function SolveInhomEVP(A::Array, B::Array, c::Array, d::Array; K₀=Nothing, a₀=Nothing,
		tol::Number=1e-6, method::Int=0, m::Int=2, sqg::Bool=false, warn::Bool=true)
	
	# Ensure that m is set correctly for SQG case

	if sqg
		
		m = 1
		
	end

	# Ensure K₀ is a Vector

	if K₀ isa Number
		
		K₀ = [K₀]
		
	end

	# Determine number of layers and number of coefficients from input
	
	N = size(d)[2]
	M = Int(size(d)[1]/N)

	# Use root finding method for multi-layer systems

	if N > 1
		
		method = 1
		
	end

	if method == 0
	
		# Use analytical solution if method = 0

		# Set K₀ value if none given
		
		if K₀ == Nothing
			
			K₀ = [4]
			
		end

		# Reformat inputs as arrays with correct shape

		B, O = reshape(B, M, M), zeros(M, M)
		dᵀ, c₀, c₁ = permutedims(d), c[:, 1], c[:, 2]

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

		v = abs.((λ .- K₀.^m).^2)
		i = argmin(v[.!isnan.(v)])

		# Define K as requested eigenvalue and invert system for a

		K = reshape([λ[i]], 1, 1).^(1/m)
		a = reshape((A - K.^m .* B) \ (c₀ + K.^m .* c₁), M, 1)

	end

	if method == 1

		# Use root finding method if method = 1

		# Calculate basis vectors spanning the space orthogonal to the d vectors
	
		e, V, iₑ = OrthogSpace(d)

		# Define a₀ if none given and reshape if given
	
		if a₀ == Nothing
			
			a₀ = vcat(-10*ones(N, 1), zeros(N*(M-1), 1))
			
		else
			
			a₀ = reshape(permutedims(a₀), N*M, 1)
			
		end

		# Define K₀ if none given and reshape if given

		if K₀ == Nothing
			
			K₀ = 5*ones(N, 1)
			
		else
			
			K₀ = reshape(K₀, N, 1)
			
		end

		# Calculate x₀ by projecting a₀ onto the new space

		x₀ = V \ a₀

		# Define x₀ as vector of K values and a₀ values in the new space

		x₀ = vcat(K₀.^m, x₀[iₑ])

		# Define function and derivative for root finding

		fj! = (F, J, x) -> InhomEVP_F!(F, J, x, A, B, c, e)

		# Solve for F(x) = 0

		x = nlsolve(only_fj!(fj!), x₀, ftol=tol).zero

		# Extract K values from first few entries of x

		K = (complex(reshape(x[1:N], 1, N))).^(1/m)

		# Calculate a values by projecting x back to original space

		a = permutedims(reshape(e * x[N+1:N*M], N, M))

	end

	# Raise warning if root finding has converged to (possibly) wrong solution

	if imag(K) != zeros(1, N)

		if warn
		
			@warn "Solution has complex K, generally corresponding passive layers."

		end
		
	end

	# Remove any imaginary parts of K and ignore coefficients less than `tol` value

	K = real(K)
	a[abs.(a) .< tol] .= 0

	return K, a

end

"""
Function: `InhomEVP_F!(F, J, x, A, B, c, d, e)`

Calculates the function F and it's derivatives, J, at a given point x

Arguments:
 - `F`, `J`: values of F and it's derivatives, updated by function
 - `x`: evaluation point, Array
 - `A`, `B`, `c`: inhomogeneous eigenvalue problem terms, Arrays
 - `e`: basis spanning the space perpendicular to the d[n], Array
"""
function InhomEVP_F!(F, J, x::Array, A::Array, B::Array, c::Array, e::Array)

	# Get problem size from inputs

	N, j = size(e)

	# Calculate a by projecting x into original space

	a = e * x[N-j+1:N]

	# Build matrix and RHS vector for F(x)

	M = A
	v = c[:, 1]
		
	for i in 1:N-j
		
		M = M - x[i] * B[:, :, i]
		v = v + x[i] * c[:, i + 1]
		
	end

	# Calculate and update Jacobian matrix of derivatives

	if !(J == nothing)

		for i in 1:N-j
			
			J[:, i] = -B[:, :, i] * a - c[:, i + 1]
			
		end

		J[:, N-j+1:end] .= M * e
		
	end

	# Calculate and update function value
	
	if !(F == nothing)
		
		F[:] .= M * a - v

	end

end

"""
Function: `OrthogSpace(v)`

Extends the input to an orthonormal basis over R^n using the Gram-Schmidt method

Arguments:
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

	# Orthogonalise full basis

	for j in 1:N
		
		for i in 1:j-1
			
			B[:, j] = B[:, j] - B[:, i] * dot(B[:, i], B[:, j]) / norm(B[:, i])^2
			
		end
		
	end

	# Normalise all vectors
	
	B = B ./ sqrt.(sum(abs2, B, dims=1))

	# Define e as the new vectors added to the basis (B is the full basis)

	e = B[:, iₑ]
	
	return e, B, iₑ

end

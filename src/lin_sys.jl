# Contains functions for building linear eigenvalue system

# Calculates the terms A, B, c, d in the linear system

function BuildLinSys(M, λ, μ; tol=1e-6)

	N = length(μ)

	A = zeros(N*M, N*M)
	B₀ = zeros(N*M, N*M)

	for j in 0:M-1
		for k in 0:M-1

			A[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> A_func(x, λ, μ), j, k, tol)[1]
			B₀[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> B_func(x, λ, μ), j, k, tol)[1]

		end
	end

	B = zeros(N*M, N*M, N)
	c = zeros(N*M, N+1)
	d = zeros(N*M, N)
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

function ApplyPassiveLayers(A, B, c, d, ActiveLayers)

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

function IncludePassiveLayers(K, a, ActiveLayers)
	
	M, N = size(a)[1], length(ActiveLayers)

	K₁, a₁ = zeros(1, N), zeros(M, N)
	
	i = BitArray{1}(ActiveLayers)

	K₁[i] .= K
	a₁[:, i] .= a
	
	return K₁, a₁

end

function SolveInhomEVP(A, B, c, d; K₀=0, a₀=0, tol=1e-6)

	e, V, iₑ = OrthogSpace(d)

	N = size(d)[2]
	M = Int(size(d)[1]/N)
	
	if a₀ == 0
		a₀ = vcat(-10*ones(N, 1), zeros(N*(M-1), 1))
	end

	if K₀ == 0
		K₀ = 5*ones(N, 1)
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

function InhomEVP_F!(F, J, x, A, B, c, e)

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

# extends v to an orthonormal basis

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
				throw("The v must be linerly independent.")
			end
			if dot(v[:, i],B[:, iₑ[j]]) > ϵ
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
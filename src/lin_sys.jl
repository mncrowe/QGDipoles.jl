# Contains functions for building linear eigenvalue system

# Calculates the terms A, B, c, d in the linear system

function BuildLinSys(M, λ, μ, tol=1e-6)

	N = length(μ)

	A = zeros(N*M, N*M)
	Bi = zeros(N*M, N*M)

	for j in 0:M-1
		for k in 0:M-1

			A[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> A_func(x, λ, μ), j, k, tol)[1]
			Bi[j*N.+(1:N), k*N.+(1:N)] .= JJ_int(x -> B_func(x, λ, μ), j, k, tol)[1]

		end
	end

	B = zeros(N*M, N*M, N)
	c = zeros(N*M, N+1)
	d = zeros(N*M, N)
	c₀ = [ones(N); zeros((M-1)*N)]
	
	for n in 1:N
		
		K = kron(I(M), diagm((1:N).==n))
		B[:, :, n] = K * Bi
		c[:, 1] = c[:, 1] + μ[n] * (K * c₀) / 4
		c[:, n+1] = (K * c₀) / 4
		d[:, n] = kron((-1).^(0:M-1), (1:N).==n)
		
	end

	return A, B, c, d

end
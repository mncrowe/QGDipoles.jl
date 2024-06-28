# Calculates the A and B matrix terms using double Bessel integrals

# Evaluates the matrix function A(ξ, λ, μ) = K(ξ) [K(ξ) + D(μ)]⁻¹ ξ⁻¹

function A_func(ξ, λ, μ)
	
	N = length(μ)
	
	if N == 1
		
		K = @. ξ^2 + λ^2
		A = @. K / (K + μ) / ξ

	else
		
		K = I(N) * ξ^2 + diagm([λ[1]^2; 2*λ[2:end-1].^2; λ[end]^2]) +
			-diagm(1 => λ[1:end-1].^2, -1 => λ[2:end].^2)
		A = (K / (K .+ diagm(μ))) / ξ

	end	

	return A

end

# Evaluates the matrix function B(ξ, λ, μ) = [K(ξ) + D(μ)]⁻¹ ξ⁻¹

function B_func(ξ, λ, μ)
	
	N = length(μ)

	if N == 1

		K = @. ξ^2 + λ^2
		B = @. 1 / (K + μ) / ξ

	else
		
		K = I(N) * ξ^2 + diagm([λ[1]^2; 2*λ[2:end-1].^2; λ[end]^2]) +
			-diagm(1 => λ[1:end-1].^2, -1 => λ[2:end].^2)
		B = inv(K .+ diagm(μ)) / ξ

	end

	return B

end

# Evalulates the integral I = int_0^inf F(x) J_{2j+2}(x) J_{2k+2}(x) dx

function JJ_int(F, j, k, tol=1e-6)
	
	#N = size(F(0))[1]	# dimension of image of F
	d = 1e3			# domain splitting parameter
	D = 100			# domain limit for exponential contribution
	atol = 1e-4*tol		# absolute error tolerance for quadgk
	
	J(n, x) = besselj(n, x)
	Y(n, x) = bessely(n, x)
	H₁(n, x) = hankelh1(n, x)
	H₂(n, x) = hankelh2(n, x)
	
	HH₁₂(m,n,x) = @. 2*(J(m, x) * J(n, x) + Y(m, x) * Y(n, x))
	ϕ₁₁(m,n,x) = @. H₁(m, x) * H₁(n, x) / exp(2*im*x)
	ϕ₂₂(m,n,x) = @. H₂(m, x) * H₂(n, x) / exp(-2*im*x)
	
	F₁(x) = F(x) .* J(2*j+2,x) * J(2*k+2,x)
	F₂(x) = F(x) .* HH₁₂(2*j+2, 2*k+2, x)/4
	F₃(x) = F(d + im*x) .* (ϕ₁₁(2*j+2, 2*k+2, d + im*x) * exp(-2*x))
	F₄(x) = F(d - im*x) .* (ϕ₂₂(2*j+2, 2*k+2, d - im*x) * exp(-2*x))
	
	I₁ = quadgk(F₁ ,0, d, rtol=tol, atol=atol)
	I₂ = quadgk(F₂, d, Inf, rtol=tol, atol=atol)
	I₃ = quadgk(F₃, 0, D, rtol=tol, atol=atol)
	I₄ = quadgk(F₄, 0, D, rtol=tol, atol=atol)
	
	I = I₁[1] + I₂[1] - imag(exp(2*im*d)*I₃[1] - exp(-2*im*d)*I₄[1])/4
	ϵ = I₁[2] + I₂[2] + I₃[2] + I₄[2]

	return I, ϵ

end
# This file contains functions for creating the modon solution for a given grid and parameters

# Define the Zernike radial function R_n using the Jacobi polynomial

function ZernikeR(n,x)
	
	y = @. (-1)^n * x * jacobi(2*x^2 - 1, n, 0, 1) * (x <= 1)

	return y
end

# grid structure for grid

struct GridStruct
	x
	y
	kr::Array{Float64}
	l::Array{Float64}
	Krsq::Array{Float64}
end

# define a grid consisting of a physical space grid and a Fourier grid for the FFTW real transform

function CreateGrid(Nx, Ny, Lx, Ly)

	if length(Lx) == 2
	
		x = range(Lx[1], step = (Lx[2] - Lx[1]) / Nx, length = Nx)
		kr = Array(reshape(rfftfreq(Nx, 2π/(Lx[2]-Lx[1])*Nx), (Int(Nx/2 + 1), 1)))

	end

	if length(Lx) == 1
	
		x = range(-Lx/2, step = Lx / Nx, length = Nx)
		kr = Array(reshape(rfftfreq(Nx, 2π/Lx*Nx), (Int(Nx/2 + 1), 1)))

	end

	if length(Ly) == 2
	
		y = range(Ly[1], step = (Ly[2] - Ly[1]) / Ny, length = Ny)
		l =  Array(reshape( fftfreq(Ny, 2π/(Ly[2]-Ly[1])*Ny), (1, Ny)))

	end

	if length(Ly) == 1
	
		y = range(-Ly/2, step = Ly / Ny, length = Ny)
		l =  Array(reshape( fftfreq(Ny, 2π/Ly*Ny), (1, Ny)))

	end

	Krsq = @. kr^2 + l^2

	return GridStruct(x, y, kr, l, Krsq)

end

# function for calculating ψ and q in each layer using coefficients and vortex parameters

function Calc_ψq(a, U, ℓ, R, β, grid, x₀=[0, 0])
	
	M, N = size(a)
	Nx, Ny = length(grid.x), length(grid.y)
	
	x, y = reshape(Array(grid.x), :, 1), reshape(Array(grid.y), 1, :)
	r, θ = @. sqrt((x-x₀[1])^2 + (y-x₀[2])^2), @. atan(y-x₀[2], x-x₀[1])

	F = zeros(Nx, Ny, N)
	ψh = Array{ComplexF64}(zeros(Int(Nx/2+1), Ny, N))
	qh = Array{ComplexF64}(zeros(Int(Nx/2+1), Ny, N))

	for j in 1:M
		for i in 1:N
			F[:, :, i] += a[j, i] * ZernikeR(j-1, r/ℓ);
		end
	end

	F = -U/ℓ * F .* sin.(θ)
	Fh = rfft(F, [1, 2])
	Fh[1, 1, :] .= 0

	ΔN_i = stack(inv, eachslice(ΔNCalc(grid.Krsq, R, β, U), dims=(3,4)))
	ΔN_0 = ΔNCalc(grid.Krsq, R, 0)

	for n in 1:N
		for j in 1:N
			ψh[:, :, n] .+= ΔN_i[n, j, :, :] .* Fh[:, :, j]
		end
	end

	for n in 1:N
		for j in 1:N
			qh[:, :, n] .+= ΔN_0[n, j, :, :] .* ψh[:, :, j]
		end
	end

	ψ = irfft(ψh, Nx, [1, 2])
	q = irfft(qh, Nx, [1, 2])
	
	return ψ, q

end

# defines the ΔN matrix used to invert for ψ and q

function ΔNCalc(K², R, β, U=1)
	
	N = length(R)
	Nk, Nl = size(K²)
	ϵ = 1e-16
	
	if length(β) < N
		βU⁻¹ = zeros(N)
	else
		βU⁻¹ = β / U
	end

	if N == 1

		ΔN = -reshape(K² .+ (1/R^2 + βU⁻¹ + ϵ), 1, 1, Nk, Nl)

	else

		ΔN = -diagm([R[1]^-2; 2*R[2:end-1].^-2; R[end]^-2]) - diagm(βU⁻¹) +
			diagm(1 => R[1:end-1].^-2, -1 => R[2:end].^-2) .-
			I(N) .* reshape(K² .+ ϵ, 1, 1, Nk, Nl)

	end

	return ΔN

end
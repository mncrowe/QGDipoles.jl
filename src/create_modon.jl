"""
This file contains the functions to calculate the streamfunction and vorticity
fields in each layer given a set of coefficients. These are determined by first
calculating:

F = -U/ℓ sin θ sum_{j=0}^{M-1} a_j R_j(r/ℓ),

where (U, ℓ) denotes the vortex speed and radius, a_j is the j^th column of the
coefficient matrix `a`, M is the number of coefficients calculated and (r, θ)
are plane polar coordinates about the vortex center. The function R_j denotes
the Zernike radial function of order j. The streamfunction and vorticity can be
calculated from F using:

Δ_N(β) ψ = F,

and

q = Δ_N(0) ψ,

where

Δ_N(β) = [∇²-R[1]⁻²-β[1]/U,      R[1]⁻²     ,   0   , ...             ...        0;
	       R[2]⁻²     , ∇²-R[1]⁻²-β[1]/U, R[2]⁻², ...             ...        0;
	        ...       ,        ...      ,  ...  , ...   ... ,     ...         ;
	         0        ,         0       ,   0   , ... R[N]⁻², ∇²-R[N]⁻²-β[N]/U]

We calculate ψ and q in Fourier space where ∇² = -(k² + l²) using the FFTW package.
Since ψ and q are real, we use the `rfft` and `irfft` functions and create wavenumber
grids which are consistent with this choice.

This package contains a function for creating the grid (`CreateGrid`) however these
functions are designed to be consistent with the `TwoDGrid` from FourierFlows.jl and
GeophysicalFlows.jl so will also work with this grid.
"""


"""
Function: ZernikeR

Define the Zernike radial function using the jacobi function from SpecialFunctions

Arguments:
 - n: order, Integer
 - x: evaluation point, Number or Array
"""

function ZernikeR(n::Int, x::Union{Number,Array})
	
	y = @. (-1)^n * x * jacobi(2*x^2 - 1, n, 0, 1) * (x <= 1)

	return y
end

"""
Structure: GridStruct

Stores the grid variables in physical and Fourier space

Arguments:
 - x, y: x and y points in physical space, Ranges
 - kr, l: x and y points in Fourier space, Arrays
 - Krsq: kr²+l² in Fourier space, Array
"""

struct GridStruct
	x
	y
	kr::Array{Float64}
	l::Array{Float64}
	Krsq::Array{Float64}
end

"""
Function: CreateGrid

Define the numerical grid as a `GridStruct`

Arguments:
 - Nx, Ny: number of gridpoints in x and y directions, Integers
 - Lx, Ly: x and y domains, either vectors of endpoints or lengths, Numbers or Vectors
"""

function CreateGrid(Nx::Int, Ny::Int, Lx::Union{Number,Vector}, Ly::Union{Number,Vector})

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

"""
Function: Calc_ψq

Calculate ψ and q in each layer using coefficients and vortex parameters

Arguments:
 - a: M x N array of coefficients, Array
 - (U, ℓ): vortex speed and radius, Numbers
 - (R, β): Rossby radii and (y) PV gradients in each layer, Numbers or Vectors
 - grid: grid structure containing x, y, and Krsq
 - x₀: position of vortex center, vector (default: [0, 0])
"""

function Calc_ψq(a::Array, U::Number, ℓ::Number, R::Union{Number,Vector}, β::Union{Number,Vector}, grid, x₀::Vector=[0, 0])
	
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


"""
Function: ΔNCalc

Defines the Δ_N(β) matrix used to invert for ψ and q

Arguments:
 - K²: value of k²+l² in Fourier space, Array
 - (R, β): Rossby radii and (y) PV gradients in each layer, Numbers or Vectors
 - U: vortex speed, Number (default: 1)
"""

function ΔNCalc(K²::Array, R::Union{Number,Vector}, β::Union{Number,Vector}, U::Number=1)
	
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
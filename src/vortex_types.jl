"""
This file contains functions which build vortices and vortex parameter sets as structures.

These functions are intended to allow all fields and variables associated with a given
solution to be included in the same structure and allow parameters to be checked.

"""

"""
Structure: `LQGParams`

Stores the parameters for an LQG dipolar vortex solution

Arguments:
 - `U`: Vortex speed
 - `ℓ`: Vortex radius
 - `R`: Rossby radius
 - `β`: background PV gradient
 - `ActiveLayers`: 1 => layer contains vortex region
 - `H`: thickness of each layer
 - `x₀`: Vortex position
 - `α`: Direction of vortex propagation
 - `M`: number of coefficients in Zernike expansion
 - `tol`: maximum error in solution evaluation
 - `K₀`: initial guess for eigenvalue
 - `a₀`: initial guess for coefficients
 - `UseAnalytic`: use analytic solution (1-layer only)
 - `CalcVelocity`: flag to determine if velocity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
 - `CalcEnstrophy`: flag to determine if enstrophy is calculated
"""
struct LQGParams
	U::Number
	ℓ::Number
	R::Union{Number,Vector}
	β::Union{Number,Vector}
	ActiveLayers::Union{Number,Vector}
	H::Union{Number,Vector}
	x₀::Vector
	α::Number
	M::Int
	tol::Number
	K₀::Union{Number,Array,Nothing}
	a₀::Union{Array,Nothing}
	UseAnalytic::Bool
	CalcVelocity::Bool
	CalcEnergy::Bool
	CalcEnstrophy::Bool
end

"""
Structure: `SQGParams`

Stores the parameters for an SQG dipolar vortex solution

Arguments:
 - `U`: Vortex speed
 - `ℓ`: Vortex radius
 - `R`: Rossby radius
 - `β`: background PV gradient
 - `x₀`: Vortex position
 - `α`: Direction of vortex propagation
 - `M`: number of coefficients in Zernike expansion
 - `tol`: maximum error in solution evaluation
 - `K₀`: initial guess for eigenvalue
 - `a₀`: initial guess for coefficients
 - `CalcVelocity`: flag to determine if velocity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
 
"""
struct SQGParams
	U::Number
	ℓ::Number
	R::Vector
	β::Number
	x₀::Vector
	α::Number
	M::Int
	tol::Number
	K₀::Union{Number,Array,Nothing}
	a₀::Union{Array,Nothing}
	CalcVelocity::Bool
	CalcEnergy::Bool
end

"""
Structure: `LQGVortex`

Stores fields and diagnostics for an LQG dipolar vortex solution

Arguments:
 - `params`: Vortex params
 - `ψ`: streamfunction
 - `q`: potential vorticity anomaly
 - `K`: eigenvalue
 - `a`: coefficient matrix
 - `u`: x velocity
 - `v`: y velocity
 - `KE`: kinetic energy
 - `PE`: potential energy
 - `EN`: enstrophy 
"""
struct LQGVortex
	params::LQGParams
	ψ::Union{CuArray,Array}
	q::Union{CuArray,Array}
	K::Array
	a::Union{Array,Nothing}
	u::Union{CuArray,Array,Nothing}
	v::Union{CuArray,Array,Nothing}
	KE::Union{Vector,Nothing}
	PE::Union{Vector,Nothing}
	EN::Union{Vector,Nothing}
end

"""
Structure: `SQGVortex`

Stores fields and diagnostics for an SQG dipolar vortex solution

Arguments:
 - `params`: Vortex params
 - `ψ`: surface streamfunction
 - `b`: surface buoyancy
 - `K`: eigenvalue
 - `a`: coefficient matrix
 - `u`: x velocity
 - `v`: y velocity
 - `E`: domain integrated energy
 - `SPE`: surface potential energy
"""
struct SQGVortex
	params::SQGParams
	ψ::Union{CuArray,Array}
	b::Union{CuArray,Array}
	K::Array
	a::Array
	u::Union{CuArray,Array,Nothing}
	v::Union{CuArray,Array,Nothing}
	E::Union{Vector, Nothing}
	SPE::Union{Vector, Nothing}
end

"""
Function: `DefLQGParams`

Defines an LQGParams structure using the given inputs

Arguments:
 - `U`: Vortex speed
 - `ℓ`: Vortex radius
 - `R`: Rossby radius
 - `β`: background PV gradient
 - `ActiveLayers`: 1 => layer contains vortex region
 - `H`: thickness of each layer
 - `x₀`: Vortex position
 - `α`: Direction of vortex propagation
 - `M`: number of coefficients in Zernike expansion
 - `tol`: maximum error in solution evaluation
 - `K₀`: initial guess for eigenvalue
 - `a₀`: initial guess for coefficients
 - `UseAnalytic`: use analytic solution (1-layer only)
 - `CalcVelocity`: flag to determine if velocity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
 - `CalcEnstrophy`: flag to determine if enstrophy is calculated
"""
function DefLQGParams(U::Number=1,
		      ℓ::Number=1,
		      R::Union{Number,Vector}=Inf,
		      β::Union{Number,Vector}=0,
	   ActiveLayers::Union{Number,Vector}=1,
		      H::Union{Number,Vector}=1,
		     x₀::Vector=[0, 0],
		      α::Number=0,
		      M::Int=8;
	            tol::Number=1e-6,
		     K₀::Union{Number,Array,Nothing}=nothing,
		     a₀::Union{Array,Nothing}=nothing,
	    UseAnalytic::Bool=false,
	   CalcVelocity::Bool=false,
             CalcEnergy::Bool=false,
	  CalcEnstrophy::Bool=false)
	
	return LQGParams(U, ℓ, R, β, ActiveLayers, H, x₀, α, M,
		tol, K₀, a₀, UseAnalytic, CalcVelocity, CalcEnergy, CalcEnstrophy)
end

"""
Function: `DefSQGParams`

Defines an SQGParams structure using the given inputs

Arguments:
 - `U`: Vortex speed
 - `ℓ`: Vortex radius
 - `R`: Rossby radius
 - `β`: background PV gradient
 - `x₀`: Vortex position
 - `α`: Direction of vortex propagation
 - `M`: number of coefficients in Zernike expansion
 - `tol`: maximum error in solution evaluation
 - `K₀`: initial guess for eigenvalue
 - `a₀`: initial guess for coefficients
 - `CalcVelocity`: flag to determine if velocity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
"""
function DefSQGParams(U::Number=1,
		      ℓ::Number=1,
		      R::Vector=[Inf, Inf],
		      β::Number=0,
		     x₀::Vector=[0, 0],
		      α::Number=0,
		      M::Int=12;
	            tol::Number=1e-6,
		     K₀::Union{Number,Array,Nothing}=nothing,
		     a₀::Union{Array,Nothing}=nothing,
	   CalcVelocity::Bool=false,
             CalcEnergy::Bool=false)

	return SQGParams(U, ℓ, R, β, x₀, α, M, tol, K₀, a₀, CalcVelocity, CalcEnergy)
end

"""
Function: `DefLQGVortex`

Defines an LQGVortex solution structure using the given inputs

Arguments:
 - `grid`: grid structure
 - `U`: Vortex speed
 - `ℓ`: Vortex radius
 - `R`: Rossby radius
 - `β`: background PV gradient
 - `ActiveLayers`: 1 => layer contains vortex region
 - `H`: thickness of each layer
 - `x₀`: Vortex position
 - `α`: Direction of vortex propagation
 - `M`: number of coefficients in Zernike expansion
 - `tol`: maximum error in solution evaluation
 - `K₀`: initial guess for eigenvalue
 - `a₀`: initial guess for coefficients
 - `UseAnalytic`: use analytic solution (1-layer only)
 - `CalcVelocity`: flag to determine if velocity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
 - `CalcEnstrophy`: flag to determine if enstrophy is calculated
"""
function DefLQGVortex(grid,
		      U::Number=1,
		      ℓ::Number=1,
		      R::Union{Number,Vector}=Inf,
		      β::Union{Number,Vector}=0,
	   ActiveLayers::Union{Number,Vector}=1,
		      H::Union{Number,Vector}=1,
		     x₀::Vector=[0, 0],
		      α::Number=0,
		      M::Int=8;
	            tol::Number=1e-6,
		     K₀::Union{Number,Array,Nothing}=nothing,
		     a₀::Union{Array,Nothing}=nothing,
	    UseAnalytic::Bool=false,
	   CalcVelocity::Bool=false,
             CalcEnergy::Bool=false,
	  CalcEnstrophy::Bool=false)

	params = DefLQGParams(U, ℓ, R, β, ActiveLayers, H, x₀, α, M;
			tol, K₀, a₀, UseAnalytic, CalcVelocity, CalcEnergy, CalcEnstrophy)

	N = length(R)

	if UseAnalytic
		if N == 1
			ψ, q, K = CreateLRD(grid, U, ℓ, R, β, x₀, α)
			a = nothing
		else
			@error "UseAnalytic = true not supported for N > 1"
		end
	else
		ψ, q, K, a = CreateModonLQG(grid, M, U, ℓ, R, β, ActiveLayers, x₀, α; K₀, a₀, tol)
	end

	if CalcVelocity
		u, v = Calc_uv(ψ, grid)
	else
		u, v = nothing, nothing
	end

	if CalcEnergy
		KE, PE = EnergyLQG(grid, ψ, R, H)
	else
		KE, PE = nothing, nothing
	end

	if CalcEnstrophy
		EN = EnstrophyLQG(grid, q, H)
	else
		EN = nothing
	end

	return LQGVortex(params, ψ, q, K, a, u, v, KE, PE, EN)
end

"""
Function: `DefSQGVortex`

Defines an SQGVortex solution structure using the given inputs

Arguments:
 - `grid`: grid structure
 - `U`: Vortex speed
 - `ℓ`: Vortex radius
 - `R`: Rossby radius
 - `β`: background PV gradient
 - `x₀`: Vortex position
 - `α`: Direction of vortex propagation
 - `M`: number of coefficients in Zernike expansion
 - `tol`: maximum error in solution evaluation
 - `K₀`: initial guess for eigenvalue
 - `a₀`: initial guess for coefficients
 - `CalcVelocity`: flag to determine if velocity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
"""
function DefSQGVortex(grid,
		      U::Number=1,
		      ℓ::Number=1,
		      R::Vector=[Inf, Inf],
		      β::Number=0,
		     x₀::Vector=[0, 0],
		      α::Number=0,
		      M::Int=12;
	            tol::Number=1e-6,
		     K₀::Union{Number,Array,Nothing}=nothing,
		     a₀::Union{Number,Array,Nothing}=nothing,
	   CalcVelocity::Bool=false,
             CalcEnergy::Bool=false)

	params = DefSQGParams(U, ℓ, R, β, x₀, α, M; tol, K₀, a₀, CalcVelocity, CalcEnergy)

	ψ, b, K, a = CreateModonSQG(grid, M, U, ℓ, R, β, x₀, α; K₀, a₀, tol)

	if CalcVelocity
		u, v = Calc_uv(ψ, grid)
	else
		u, v = nothing, nothing
	end

	if CalcEnergy
		E, SPE = EnergySQG(grid, ψ, b, R[2])
	else
		E, SPE = nothing, nothing
	end

	return SQGVortex(params, ψ, b, K, a, u, v, E, SPE)
end

"""
Function: `DefLQGVortex`

Defines an LQGVortex solution structure using the given inputs

Arguments:
 - `grid`: grid structure
 - `params`: vortex parameters, LQGParams structure
"""
function DefLQGVortex(grid, params::LQGParams)

	return DefLQGVortex(grid, params.U, params.ℓ, params.R, params.β, params.ActiveLayers, params.H,
			params.x₀, params.α, params.M; params.tol, params.K₀, params.a₀, params.UseAnalytic,
			params.CalcVelocity, params.CalcEnergy, params.CalcEnstrophy)
end

"""
Function: `DefSQGVortex`

Defines an SQGVortex solution structure using the given inputs

Arguments:
 - `grid`: grid structure
 - `params`: vortex parameters, LQGParams structure
"""
function DefSQGVortex(grid, params::SQGParams)

	return DefSQGVortex(grid, params.U, params.ℓ, params.R, params.β, params.x₀, params.α, params.M;
			params.tol, params.K₀, params.a₀, params.CalcVelocity, params.CalcEnergy)
end

"""
Base.summary function for custom type `LQGParams`
"""
function Base.summary(g::LQGParams)
	return string("Parameter set structure for an LQG vortex solution (LQGParams)")
end


"""
Base.summary function for custom type `SQGParams`
"""
function Base.summary(g::SQGParams)
	return string("Parameter set structure for an SQG vortex solution (SQGParams)")
end


"""
Base.summary function for custom type `LQGVortex`
"""
function Base.summary(g::LQGVortex)
	return string("Vortex solution structure for an LQG model (LQGVortex)")
end


"""
Base.summary function for custom type `SQGVortex`
"""
function Base.summary(g::SQGVortex)
	return string("Vortex solution structure for an SQG model (SQGVortex)")
end

"""
Base.show function for custom type `LQGParams`
"""
function Base.show(io::IO, p::LQGParams)

	N = length(p.R)
	
	a₀_given = ~(p.a₀ isa Nothing)

	return print(io, "LQGParams\n",
               "  ├───────── number of layers (N): ", N, "\n",
               "  ├───────────── vortex speed (U): ", p.U, "\n",
               "  ├──────────── vortex radius (ℓ): ", p.ℓ, "\n",
               "  ├──────────── Rossby radius (R): ", p.R, "\n",
               "  ├───────── PV gradient in y (β): ", p.β, "\n",
               "  ├─ active layers (ActiveLayers): ", p.ActiveLayers, "\n",
               "  ├───────────── layer depths (H): ", p.H, "\n",
               "  ├───────── vortex position (x₀): ", p.x₀, "\n",
               "  ├───────────── vortex angle (α): ", p.α, "\n",
               "  ├─── number of coefficients (M): ", p.M, "\n",
               "  ├──────── error tolerance (tol): ", p.tol, "\n",
               "  ├──── guess for eigenvalue (K₀): ", p.K₀, "\n",
               "  ├── guess for coeffs (a₀) given: ", a₀_given, "\n",
               "  ├──────── use analytic solution: ", p.UseAnalytic, "\n",
               "  ├─────────── calculate velocity: ", p.CalcVelocity, "\n",
               "  ├───────────── calculate energy: ", p.CalcEnergy, "\n",
               "  └────────── calculate enstrophy: ", p.CalcEnstrophy)    

end

"""
Base.show function for custom type `SQGParams`
"""
function Base.show(io::IO, p::SQGParams)
	
	a₀_given = ~(p.a₀ isa Nothing)

	return print(io, "SQGParams\n",
               "  ├───────────── vortex speed (U): ", p.U, "\n",
               "  ├──────────── vortex radius (ℓ): ", p.ℓ, "\n",
               "  ├──────────── Rossby radius (R): ", p.R, "\n",
               "  ├───────── PV gradient in y (β): ", p.β, "\n",
               "  ├───────── vortex position (x₀): ", p.x₀, "\n",
               "  ├───────────── vortex angle (α): ", p.α, "\n",
               "  ├─── number of coefficients (M): ", p.M, "\n",
               "  ├──────── error tolerance (tol): ", p.tol, "\n",
               "  ├──── guess for eigenvalue (K₀): ", p.K₀, "\n",
               "  ├── guess for coeffs (a₀) given: ", a₀_given, "\n",
               "  ├─────────── calculate velocity: ", p.CalcVelocity, "\n",
               "  └───────────── calculate energy: ", p.CalcEnergy)    

end

"""
Base.show function for custom type `LQGVortex`
"""
function Base.show(io::IO, p::LQGVortex)

	contains_a = ~(p.a isa Nothing)
	contains_vel = ~(p.u isa Nothing)
	contains_ener = ~(p.KE isa Nothing)
	contains_enst = ~(p.EN isa Nothing)

	if p.ψ isa CuArray
		dev = "GPU"
	else
		dev = "CPU"
	end

	UseAnalytic = p.params.UseAnalytic

	return print(io, "LQGVortex\n",
               "  ├─ params structure (LQGParams): varname.params\n",
               "  ├─────────────────────── device: ", dev, "\n",
               "  ├──────── use analytic solution: ", UseAnalytic, "\n",
               "  ├────────── contains ψ, q and K: ", true, "\n",
               "  ├─────────────────── contains a: ", contains_a, "\n",
               "  ├────────── contains velocities: ", contains_vel, "\n",
               "  ├──────────── contains energies: ", contains_ener, "\n",
               "  └─────────── contains enstrophy: ", contains_enst)    

end

"""
Base.show function for custom type `SQGVortex`
"""
function Base.show(io::IO, p::SQGVortex)

	contains_vel = ~(p.u isa Nothing)
	contains_ener = ~(p.E isa Nothing)

	if p.ψ isa CuArray
		dev = "GPU"
	else
		dev = "CPU"
	end

	return print(io, "SQGVortex\n",
               "  ├─ params structure (SQGParams): varname.params\n",
               "  ├─────────────────────── device: ", dev, "\n",
               "  ├─────── contains ψ, b, K and a: ", true, "\n",
               "  ├────────── contains velocities: ", contains_vel, "\n",
               "  └──────────── contains energies: ", contains_ener)    

end
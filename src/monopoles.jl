"""
While primarily intended for creating dipolar vortex solutions, this package contains
additional functions which can be used to generate some simple monopolar vortex solutions
too. These functions are included in this file.
"""

"""
Function: `CreateRankine(grid, ℓ=1, x₀=[0, 0])`

Calculates the Rankine vortex for a 1 layer system

Arguments:
 - `grid`: grid structure containing x, y, and Krsq
 - `ℓ`: vortex speed and radius, Numbers (default: `1`)
 - `Γ`: vortex circulation (default: `2π`)
 - `x₀`: position of vortex center, vector (default: `[0, 0]`)

Note: This function outputs (u, v) directly since the solution has discontinuous velocity at
the vortex boundary, r = ℓ.
"""
function CreateRankine(grid, ℓ::Number=1, Γ::Number=2π, x₀::Vector=[0, 0])

	# Create Cartesian and polar grids

	x, y = reshape(Array(grid.x), :, 1), reshape(Array(grid.y), 1, :)
	r, θ = @. sqrt((x-x₀[1])^2 + (y-x₀[2])^2), @. atan(y-x₀[2], x-x₀[1])

	# Calculate ψ and q using analytic result

	ψ = @.   Γ * (log(r / ℓ) * (r >= ℓ) + 1/2 * (r^2 / ℓ^2 - 1) * (r < ℓ))
	q = @.   Γ * (2 / ℓ^2) * (r < ℓ)
	u = @. - Γ * (1 / r * (r >= ℓ) + r / ℓ^2  * (r < ℓ)) * sin(θ)
	v = @.   Γ * (1 / r * (r >= ℓ) + r / ℓ^2  * (r < ℓ)) * cos(θ)

	# Move result to GPU if `cuda=true` in grid

	if grid.Krsq isa CuArray
		
		ψ, q = CuArray(ψ), CuArray(q)
		u, v = CuArray(u), CuArray(v)
		
	end
	
	return ψ, q, u, v

end




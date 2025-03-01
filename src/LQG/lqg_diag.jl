"""
This file contains functions to calculate the energy and enstrophy for dipolar vortex
solutions in the LQG model.

For the (multi-layer) LQG model we have:

``KE = \\frac{H_i}{2H} \\int_A |\\nabla\\psi_i|^2 dx dy,``

``PE = \\frac{H_i}{2H R_i^2} \\int_A |\\psi_i - \\psi_{i+1}|^2 dx dy,``

For 1 layer (equivalent barotropic) QG, the PE is given by

``PE = \\frac{1}{2 R^2} \\int_A |\\psi|^2 dx dy,``

and the KE is the same as the multi-layer case.

"""

"""
    EnergyLQG(grid, ψ, R, H=[1])

Calculates the kinetic and potential energy for the LQG system

# Arguments:
 - `grid`: grid structure containing `Krsq`
 - `ψ`: streamfunction in each layer, Array or CuArray
 - `R`: Rossby radius in each layer, Number or Vector
 - `H`: Thickness of each layer, Number or Vector
"""
function EnergyLQG(
    grid,
    ψ::Union{CuArray,Array},
    R::Union{Number,Vector},
    H::Union{Number,Vector} = [1],
)

    Nx, Ny = size(ψ)[1:2]
    N = Int(length(ψ) / (Nx * Ny))

    if N == 1

        ψh = rfft(ψ, [1, 2])
        eh = sqrt.(grid.Krsq) .* ψh

        KE = [AreaInteg2(eh, grid) / 2]
        PE = [AreaInteg2(ψh, grid) ./ R .^ 2 / 2]

    else

        D = sum(H)
        KE, PE = zeros(N), zeros(N - 1)

        ψh = rfft(ψ, [1, 2])
        eh = sqrt.(grid.Krsq) .* ψh

        for i = 1:N
            KE[i] = H[i] / (2 * D) * AreaInteg2(eh[:, :, i], grid)
        end

        for i = 1:N-1
            PE[i] = H[i] / (2 * D * R[i]^2) * AreaInteg2(ψh[:, :, i] - ψh[:, :, i+1], grid)
        end

    end

    return KE, PE

end

"""
    EnstrophyLQG(grid, q, H=[1])

Calculates the enstrophy for the LQG system

# Arguments:
 - `grid`: grid structure containing `Krsq`
 - `q`: potential vorticity anomaly in each layer, Array or CuArray
 - `H`: Thickness of each layer, Number or Vector
"""
function EnstrophyLQG(grid, q::Union{CuArray,Array}, H::Union{Number,Vector} = [1])

    Nx, Ny = size(q)[1:2]
    N = Int(length(q) / (Nx * Ny))

    D = sum(H)

    EN = zeros(N)

    for i = 1:N
        EN[i] = H[i] / (2 * D) * AreaInteg2(q[:, :, i], grid)
    end

    return EN

end
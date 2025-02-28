"""
This file contains functions to calculate the energy and enstrophy for dipolar vortex
solutions.

For the (multi-layer) LQG model we have:

``KE = \\frac{H_i}{2H} \\int_A |\\nabla\\psi_i|^2 dx dy,``

``PE = \\frac{H_i}{2H R_i^2} \\int_A |\\psi_i - \\psi_{i+1}|^2 dx dy,``

For 1 layer (equivalent barotropic) QG, the PE is given by

``PE = \\frac{1}{2 R^2} \\int_A |\\psi|^2 dx dy,``

and the KE is the same as the multi-layer case.

For SQG, the two quantities of interest are the total domain averaged energy

``E = -\\frac{1}{2} \\int_A \\psi b dx dy,``

and the surface potential energy

``SPE = \\frac{1}{2} \\int_A |b + \\psi / R^\\prime|^2 dx dy.``

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

"""
    EnergySQG(grid, ψ, b, R′)

Calculates the energies for the SQG system; the total domain integrated energy
and the surface potential energy

# Arguments:
 - `grid`: grid structure containing `Krsq`
 - `ψ`: surface streamfunction, Array or CuArray
 - `b`: surface buoyancy, , Array or CuArray
 - `R′`: reduced barotropic Rossby radius, Number (default: `Inf`)

Note: the surface potential energy is sometimes referred to as the generalised
enstrophy or the buoyancy variance.
"""
function EnergySQG(grid, ψ::Union{CuArray,Array}, b::Union{CuArray,Array}, R′::Number = Inf)

    ψh = rfft(ψ)
    bh = rfft(b)

    eh = sqrt.(ψh .* bh)
    sh = bh + ψh / R′

    E = [AreaInteg2(eh, grid) / 2]
    SPE = [AreaInteg2(sh, grid) / 2]

    return E, SPE

end

"""
    AreaInteg2(f, grid)

Calculates the integral ``I = ∫_A f^2 \\mathrm{d}A`` where ``A``
is the 2D domain described by `grid`.

# Arguments:
 - `f`: input Array in real or Fourier space
 - `grid`: grid structure

Note: f can be entered in real space or Fourier space, we use the rfft function
to calculate the Fourier transform so array sizes can distinguish the two.
"""
function AreaInteg2(f::Union{CuArray,Array}, grid::GridStruct)

    # Get grid parameters

    Nkr = length(grid.kr)
    Nx, Ny = length(grid.x), length(grid.y)
    Δx, Δy = grid.x[2] - grid.x[1], grid.y[2] - grid.y[1]

    # If input array is in real space, apply Fourier transform

    if size(f)[1:2] == (Nx, Ny)

        fh = rfft(f, [1, 2])

    else

        fh = f

    end

    # Add up components in Fourier space, non-edge elements are counted twice for an rfft

    I =
        sum(abs2, fh[1, :]) +             # kr = 0 (edge)
        sum(abs2, fh[Nkr, :]) +         # kr = Nx/2 (edge)
        2 * sum(abs2, fh[2:Nkr-1, :])  # sum twice for non-edge modes (as using rfft)

    I = I * (Δx * Δy) / (Nx * Ny) # normalization factor for fft

    return I

end

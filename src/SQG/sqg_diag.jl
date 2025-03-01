"""
This file cotains energy diagnostics for the SQG model.

For SQG, the two quantities of interest are the total domain averaged energy

``E = -\\frac{1}{2} ʃ_A ψ b dx dy,``

and the surface potential energy

``SPE = \\frac{1}{2} ʃ_A |b + ψ / R′|^2 dx dy.``

"""

"""
    EnergySQG(grid, ψ, b; R′)

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
function EnergySQG(
    grid,
    ψ::Union{CuArray,Array},
    b::Union{CuArray,Array};
    R′::Number = Inf,
)

    ψh = rfft(ψ)
    bh = rfft(b)

    eh = sqrt.(ψh .* bh)
    sh = bh + ψh / R′

    E = [AreaInteg2(eh, grid) / 2]
    SPE = [AreaInteg2(sh, grid) / 2]

    return E, SPE

end

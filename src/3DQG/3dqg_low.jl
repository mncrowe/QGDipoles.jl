"""
This file contains low-level functions for building 3D QG vortices.

...

"""




function Calc_ψq_3D(
    grid,
    a::Array;
    U::Number = 1,
    ℓ::Number = 1,
    H::Number = 1,
    N::Function = x -> 1,
    f::Number = 1,
    β::Number = 0,
    g::Number = Inf,
    s::Number = 0,
    x₀::Vector = [0, 0],
    α::Number = 0,
)

    λ = ζ -> a * f / (N(ζ * H) * H)
    μ = β * a^2 / U
    λ₀ = N(0) * sqrt(H / g)
    γ = s * H * N(-H)^2 / (f * U)

    




    ψ, q = 0, 0



    return ψ, q

end



function ΔNCalc(
    grid,


)

    ΔN = 0




    return ΔN


end
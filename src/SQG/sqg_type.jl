"""
This file contains functions which build vortices and vortex parameter sets as structures.

These functions are intended to allow all fields and variables associated with a given
solution to be included in the same structure and allow parameters to be checked.

"""

"""
    SQGParams

Stores the parameters for an SQG dipolar vortex solution

# Fields:
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
 - `CalcVorticity`: flag to determine if vorticity is calculated
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
    CalcVorticity::Bool
    CalcEnergy::Bool
end

"""
    SQGVortex

Stores fields and diagnostics for an SQG dipolar vortex solution

# Fields:
 - `params`: Vortex params
 - `ψ`: surface streamfunction
 - `b`: surface buoyancy
 - `K`: eigenvalue
 - `a`: coefficient matrix
 - `u`: x velocity
 - `v`: y velocity
 - `ζ`: vertical vorticity
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
    ζ::Union{CuArray,Array,Nothing}
    E::Union{Vector,Nothing}
    SPE::Union{Vector,Nothing}
end

"""
    DefSQGParams(; U=1, ℓ=1, R=[Inf,Inf], β=0, x₀=[0,0], α=0, M=12, tol=1e-6, K₀=nothing, a₀=nothing, CalcVelocity=false, CalcVorticity=false, CalcEnergy=false)

Defines an `SQGParams` structure using the given inputs

# Keyword arguments:
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
 - `CalcVorticity`: flag to determine if vorticity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
"""
function DefSQGParams(;
    U::Number = 1,
    ℓ::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
    x₀::Vector = [0, 0],
    α::Number = 0,
    M::Int = 12,
    tol::Number = 1e-6,
    K₀::Union{Number,Array,Nothing} = nothing,
    a₀::Union{Array,Nothing} = nothing,
    CalcVelocity::Bool = false,
    CalcVorticity::Bool = false,
    CalcEnergy::Bool = false,
)

    return SQGParams(
        U,
        ℓ,
        R,
        β,
        x₀,
        α,
        M,
        tol,
        K₀,
        a₀,
        CalcVelocity,
        CalcVorticity,
        CalcEnergy,
    )

end

"""
    UpdateParams(params::SQGParams; kwargs...)

Creates an `SQGParams` structure by replacing parameters in `params` with the given keywords

# Arguments:
 - `params`: `SQGParams` parameter structure

# Keyword arguments:
 - `kwargs...`: keyword arguments for `DefSQGParams`
"""
function UpdateParams(
    params::SQGParams;
    U::Number = params.U,
    ℓ::Number = params.ℓ,
    R::Vector = params.R,
    β::Number = params.β,
    x₀::Vector = params.x₀,
    α::Number = params.α,
    M::Int = params.M,
    tol::Number = params.tol,
    K₀::Union{Number,Array,Nothing} = params.K₀,
    a₀::Union{Array,Nothing} = params.a₀,
    CalcVelocity::Bool = params.CalcVelocity,
    CalcVorticity::Bool = params.CalcVorticity,
    CalcEnergy::Bool = params.CalcEnergy,
)

    return DefSQGParams(;
        U,
        ℓ,
        R,
        β,
        x₀,
        α,
        M,
        tol,
        K₀,
        a₀,
        CalcVelocity,
        CalcVorticity,
        CalcEnergy,
    )

end

"""
    DefSQGVortex(grid; U=1, ℓ=1, R=[Inf,Inf], β=0, x₀=[0,0], α=0, M=12, tol=1e-6, K₀=nothing, a₀=nothing, CalcVelocity=false, CalcVorticity=false, CalcEnergy=false)

Defines an `SQGVortex` solution structure using the given inputs

# Arguments:
 - `grid`: grid structure

# Keyword arguments:
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
 - `CalcVorticity`: flag to determine if vorticity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
"""
function DefSQGVortex(
    grid;
    U::Number = 1,
    ℓ::Number = 1,
    R::Vector = [Inf, Inf],
    β::Number = 0,
    x₀::Vector = [0, 0],
    α::Number = 0,
    M::Int = 12,
    tol::Number = 1e-6,
    K₀::Union{Number,Array,Nothing} = nothing,
    a₀::Union{Number,Array,Nothing} = nothing,
    CalcVelocity::Bool = false,
    CalcVorticity::Bool = false,
    CalcEnergy::Bool = false,
)

    params = DefSQGParams(;
        U,
        ℓ,
        R,
        β,
        x₀,
        α,
        M,
        tol,
        K₀,
        a₀,
        CalcVelocity,
        CalcVorticity,
        CalcEnergy,
    )

    ψ, b, K, a = CreateModonSQG(grid; U, ℓ, R, β, x₀, α, M, tol, K₀, a₀)

    if CalcVelocity
        u, v = Calc_uv(grid, ψ)
    else
        u, v = nothing, nothing
    end

    if CalcVorticity
        ζ = Calc_ζ(grid, ψ)
    else
        ζ = nothing
    end

    if CalcEnergy
        E, SPE = EnergySQG(grid, ψ, b; R′ = R[2])
    else
        E, SPE = nothing, nothing
    end

    return SQGVortex(params, ψ, b, K, a, u, v, ζ, E, SPE)
end

"""
    DefSQGVortex(grid, params)

Defines an `SQGVortex` solution structure using the given inputs

# Arguments:
 - `grid`: grid structure
 - `params`: vortex parameters, `LQGParams` structure
"""
function DefSQGVortex(grid, params::SQGParams)

    return DefSQGVortex(
        grid;
        params.U,
        params.ℓ,
        params.R,
        params.β,
        params.x₀,
        params.α,
        params.M,
        params.tol,
        params.K₀,
        params.a₀,
        params.CalcVelocity,
        params.CalcVorticity,
        params.CalcEnergy,
    )
end

"""
    UpdateVortex(grid, vortex::SQGVortex; kwargs...)

Creates an `SQGVortex` structure by replacing parameters in `vortex.params` with the given keywords

# Arguments:
 - `grid`: grid structure
 - `vortex`: `SQGVortex` structure

# Keyword arguments:
 - `kwargs...`: keyword arguments for `DefSQGParams`
"""
function UpdateVortex(grid, vortex::SQGVortex; kwargs...)

    params = UpdateParams(vortex.params; kwargs...)
    return DefSQGVortex(grid, params)

end

"""
    Base.summary

Summary method for custom type `SQGParams`
"""
function Base.summary(::SQGParams)
    return string("Parameter set structure for an SQG vortex solution (SQGParams)")
end

"""
    Base.summary

Summary method for custom type `SQGVortex`
"""
function Base.summary(::SQGVortex)
    return string("Vortex solution structure for an SQG model (SQGVortex)")
end

"""
    Base.show

Show method for custom type `SQGParams`
"""
function Base.show(io::IO, p::SQGParams)

    a₀_given = ~(p.a₀ isa Nothing)

    return print(
        io,
        "SQGParams\n",
        "  ├───────────── vortex speed (U): ",
        p.U,
        "\n",
        "  ├──────────── vortex radius (ℓ): ",
        p.ℓ,
        "\n",
        "  ├──────────── Rossby radius (R): ",
        p.R,
        "\n",
        "  ├───────── PV gradient in y (β): ",
        p.β,
        "\n",
        "  ├───────── vortex position (x₀): ",
        p.x₀,
        "\n",
        "  ├───────────── vortex angle (α): ",
        p.α,
        "\n",
        "  ├─── number of coefficients (M): ",
        p.M,
        "\n",
        "  ├──────── error tolerance (tol): ",
        p.tol,
        "\n",
        "  ├──── guess for eigenvalue (K₀): ",
        p.K₀,
        "\n",
        "  ├── guess for coeffs (a₀) given: ",
        a₀_given,
        "\n",
        "  ├─────────── calculate velocity: ",
        p.CalcVelocity,
        "\n",
        "  ├────────── calculate vorticity: ",
        p.CalcVorticity,
        "\n",
        "  └───────────── calculate energy: ",
        p.CalcEnergy,
    )

end

"""
    Base.show

Show method for custom type `SQGVortex`
"""
function Base.show(io::IO, p::SQGVortex)

    contains_vel = ~(p.u isa Nothing)
    contains_vor = ~(p.ζ isa Nothing)
    contains_ener = ~(p.E isa Nothing)

    if p.ψ isa CuArray
        dev = "GPU"
    else
        dev = "CPU"
    end

    return print(
        io,
        "SQGVortex\n",
        "  ├─ params structure (SQGParams): varname.params\n",
        "  ├─────────────────────── device: ",
        dev,
        "\n",
        "  ├─────── contains ψ, b, K and a: ",
        true,
        "\n",
        "  ├────────── contains velocities: ",
        contains_vel,
        "\n",
        "  ├─────────── contains vorticity: ",
        contains_vor,
        "\n",
        "  └──────────── contains energies: ",
        contains_ener,
    )

end

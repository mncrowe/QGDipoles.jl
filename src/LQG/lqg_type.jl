"""
This file contains functions which build vortices and vortex parameter sets as structures.

These functions are intended to allow all fields and variables associated with a given
solution to be included in the same structure and allow parameters to be checked.

"""

"""
    LQGParams

Stores the parameters for an LQG dipolar vortex solution

# Fields:
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
 - `CalcVorticity`: flag to determine if vorticity is calculated
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
    CalcVorticity::Bool
    CalcEnergy::Bool
    CalcEnstrophy::Bool
end

"""
    LQGVortex

Stores fields and diagnostics for an LQG dipolar vortex solution

# Fields:
 - `params`: Vortex params
 - `ψ`: streamfunction
 - `q`: potential vorticity anomaly
 - `K`: eigenvalue
 - `a`: coefficient matrix
 - `u`: x velocity
 - `v`: y velocity
 - `ζ`: vertical vorticity
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
    ζ::Union{CuArray,Array,Nothing}
    KE::Union{Vector,Nothing}
    PE::Union{Vector,Nothing}
    EN::Union{Vector,Nothing}
end

"""
    DefLQGParams(; U=1, ℓ=1, R=Inf, β=0, ActiveLayers=1, H=1, x₀=[0,0], α=0, M=8, tol=1e-6, K₀=nothing, a₀=nothing, UseAnalytic=false, CalcVelocity=false, CalcVorticity=false, CalcEnergy=false, CalcEnstrophy=false)

Defines an `LQGParams` structure using the given inputs

# Keyword arguments:
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
 - `CalcVorticity`: flag to determine if vorticity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
 - `CalcEnstrophy`: flag to determine if enstrophy is calculated
"""
function DefLQGParams(;
    U::Number = 1,
    ℓ::Number = 1,
    R::Union{Number,Vector} = Inf,
    β::Union{Number,Vector} = 0,
    ActiveLayers::Union{Number,Vector} = 1,
    H::Union{Number,Vector} = 1,
    x₀::Vector = [0, 0],
    α::Number = 0,
    M::Int = 8,
    tol::Number = 1e-6,
    K₀::Union{Number,Array,Nothing} = nothing,
    a₀::Union{Array,Nothing} = nothing,
    UseAnalytic::Bool = false,
    CalcVelocity::Bool = false,
    CalcVorticity::Bool = false,
    CalcEnergy::Bool = false,
    CalcEnstrophy::Bool = false,
)

    return LQGParams(
        U,
        ℓ,
        R,
        β,
        ActiveLayers,
        H,
        x₀,
        α,
        M,
        tol,
        K₀,
        a₀,
        UseAnalytic,
        CalcVelocity,
        CalcVorticity,
        CalcEnergy,
        CalcEnstrophy,
    )

end

"""
    UpdateParams(params::LQGParams; kwargs...)

Creates an `LQGParams` structure by replacing parameters in `params` with the given keywords

# Arguments:
 - `params`: `LQGParams` parameter structure

# Keyword arguments:
 - `kwargs...`: keyword arguments for `DefLQGParams`
"""
function UpdateParams(
    params::LQGParams;
    U::Number = params.U,
    ℓ::Number = params.ℓ,
    R::Union{Number,Vector} = params.R,
    β::Union{Number,Vector} = params.β,
    ActiveLayers::Union{Number,Vector} = params.ActiveLayers,
    H::Union{Number,Vector} = params.H,
    x₀::Vector = params.x₀,
    α::Number = params.α,
    M::Int = params.M,
    tol::Number = params.tol,
    K₀::Union{Number,Array,Nothing} = params.K₀,
    a₀::Union{Array,Nothing} = params.a₀,
    UseAnalytic::Bool = params.UseAnalytic,
    CalcVelocity::Bool = params.CalcVelocity,
    CalcVorticity::Bool = params.CalcVorticity,
    CalcEnergy::Bool = params.CalcEnergy,
    CalcEnstrophy::Bool = params.CalcEnstrophy,
)

    return params = DefLQGParams(;
        U,
        ℓ,
        R,
        β,
        ActiveLayers,
        H,
        x₀,
        α,
        M,
        tol,
        K₀,
        a₀,
        UseAnalytic,
        CalcVelocity,
        CalcVorticity,
        CalcEnergy,
        CalcEnstrophy,
    )

end

"""
    DefLQGVortex(grid; U=1, ℓ=1, R=Inf, β=0, ActiveLayers=1, H=1, x₀=[0,0], α=0, M=8, tol=1e-6, K₀=nothing, a₀=nothing, UseAnalytic=false, CalcVelocity=false, CalcVorticity=false, CalcEnergy=false, CalcEnstrophy=false)


Defines an `LQGVortex` solution structure using the given inputs

# Arguments:
 - `grid`: grid structure

# Keyword arguments:
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
 - `CalcVorticity`: flag to determine if vorticity is calculated
 - `CalcEnergy`: flag to determine if energy is calculated
 - `CalcEnstrophy`: flag to determine if enstrophy is calculated
"""
function DefLQGVortex(
    grid;
    U::Number = 1,
    ℓ::Number = 1,
    R::Union{Number,Vector} = Inf,
    β::Union{Number,Vector} = 0,
    ActiveLayers::Union{Number,Vector} = 1,
    H::Union{Number,Vector} = 1,
    x₀::Vector = [0, 0],
    α::Number = 0,
    M::Int = 8,
    tol::Number = 1e-6,
    K₀::Union{Number,Array,Nothing} = nothing,
    a₀::Union{Array,Nothing} = nothing,
    UseAnalytic::Bool = false,
    CalcVelocity::Bool = false,
    CalcVorticity::Bool = false,
    CalcEnergy::Bool = false,
    CalcEnstrophy::Bool = false,
)

    params = DefLQGParams(;
        U,
        ℓ,
        R,
        β,
        ActiveLayers,
        H,
        x₀,
        α,
        M,
        tol,
        K₀,
        a₀,
        UseAnalytic,
        CalcVelocity,
        CalcVorticity,
        CalcEnergy,
        CalcEnstrophy,
    )

    N = length(R)

    if UseAnalytic
        if N == 1
            ψ, q, K = CreateLRD(grid; U, ℓ, R, β, x₀, α)
            a = nothing
        else
            @error "UseAnalytic = true not supported for N > 1"
        end
    else
        ψ, q, K, a = CreateModonLQG(grid; U, ℓ, R, β, ActiveLayers, x₀, α, M, tol, K₀, a₀)
    end

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
        KE, PE = EnergyLQG(grid, ψ; R, H)
    else
        KE, PE = nothing, nothing
    end

    if CalcEnstrophy
        EN = EnstrophyLQG(grid, q; H)
    else
        EN = nothing
    end

    return LQGVortex(params, ψ, q, K, a, u, v, ζ, KE, PE, EN)
end

"""
    DefLQGVortex(grid, params)

Defines an `LQGVortex` solution structure using the given inputs

# Arguments:
 - `grid`: grid structure
 - `params`: vortex parameters, LQGParams structure
"""
function DefLQGVortex(grid, params::LQGParams)

    return DefLQGVortex(
        grid;
        params.U,
        params.ℓ,
        params.R,
        params.β,
        params.ActiveLayers,
        params.H,
        params.x₀,
        params.α,
        params.M,
        params.tol,
        params.K₀,
        params.a₀,
        params.UseAnalytic,
        params.CalcVelocity,
        params.CalcVorticity,
        params.CalcEnergy,
        params.CalcEnstrophy,
    )
end

"""
    UpdateVortex(grid, vortex::LQGVortex; kwargs...)

Creates an `LQGVortex` structure by replacing parameters in `vortex.params` with the given keywords

# Arguments:
 - `grid`: grid structure
 - `vortex`: `LQGVortex` structure

# Keyword arguments:
 - `kwargs...`: keyword arguments for `DefLQGParams`
"""
function UpdateVortex(grid, vortex::LQGVortex; kwargs...)

    params = UpdateParams(vortex.params; kwargs...)
    return DefLQGVortex(grid, params)

end

"""
    Base.summary

Summary method for custom type `LQGParams`
"""
function Base.summary(::LQGParams)
    return string("Parameter set structure for an LQG vortex solution (LQGParams)")
end

"""
    Base.summary

Summary method for custom type `LQGVortex`
"""
function Base.summary(::LQGVortex)
    return string("Vortex solution structure for an LQG model (LQGVortex)")
end

"""
    Base.show

Show method for custom type `LQGParams`
"""
function Base.show(io::IO, p::LQGParams)

    N = length(p.R)

    a₀_given = ~(p.a₀ isa Nothing)

    return print(
        io,
        "LQGParams\n",
        "  ├───────── number of layers (N): ",
        N,
        "\n",
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
        "  ├─ active layers (ActiveLayers): ",
        p.ActiveLayers,
        "\n",
        "  ├───────────── layer depths (H): ",
        p.H,
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
        "  ├──────── use analytic solution: ",
        p.UseAnalytic,
        "\n",
        "  ├─────────── calculate velocity: ",
        p.CalcVelocity,
        "\n",
        "  ├────────── calculate vorticity: ",
        p.CalcVorticity,
        "\n",
        "  ├───────────── calculate energy: ",
        p.CalcEnergy,
        "\n",
        "  └────────── calculate enstrophy: ",
        p.CalcEnstrophy,
    )

end

"""
    Base.show

Show method for custom type `LQGVortex`
"""
function Base.show(io::IO, p::LQGVortex)

    contains_a = ~(p.a isa Nothing)
    contains_vel = ~(p.u isa Nothing)
    contains_vor = ~(p.ζ isa Nothing)
    contains_ener = ~(p.KE isa Nothing)
    contains_enst = ~(p.EN isa Nothing)

    if p.ψ isa CuArray
        dev = "GPU"
    else
        dev = "CPU"
    end

    UseAnalytic = p.params.UseAnalytic

    return print(
        io,
        "LQGVortex\n",
        "  ├─ params structure (LQGParams): varname.params\n",
        "  ├─────────────────────── device: ",
        dev,
        "\n",
        "  ├──────── use analytic solution: ",
        UseAnalytic,
        "\n",
        "  ├────────── contains ψ, q and K: ",
        true,
        "\n",
        "  ├─────────────────── contains a: ",
        contains_a,
        "\n",
        "  ├────────── contains velocities: ",
        contains_vel,
        "\n",
        "  ├─────────── contains vorticity: ",
        contains_vor,
        "\n",
        "  ├──────────── contains energies: ",
        contains_ener,
        "\n",
        "  └─────────── contains enstrophy: ",
        contains_enst,
    )

end

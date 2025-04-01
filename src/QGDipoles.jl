"""
Package for creating steady modon solutions to the layered quasi-geostrophic equations and the
surface quasi geostrophic equations.

See `examples/` for example scripts and [here](https://mncrowe.github.io/QGDipoles.jl/) for documentation.
"""
module QGDipoles

# Load required packages

using Jacobi, FFTW, LinearAlgebra, SpecialFunctions, QuadGK, NLsolve, CUDA, RecipesBase

# Export names for all functions and types

export

    ## core

    SolveInhomEVP,
    #InhomEVP_F!,
    #OrthogSpace,

    #GridStruct,
    CreateGrid,
    #CartesianGrid,
    #PolarGrid,

    #JJ_int,
    #AreaInteg2,

    #ZernikeR,

    ## extras

    CreateRankine,
    Create1LMonopole,
    InvertVorticity1LQG,
    CreateLQGMonopole,
    InvertVorticityLQG,

    ## QG shared

    UpdateParams,
    UpdateVortex,

    ## LQG

    #A_func,
    #B_func,
    BuildLinSysLQG,
    ApplyPassiveLayers,
    IncludePassiveLayers,
    Calc_ψq,
    #ΔNCalc,

    CreateModonLQG,
    CreateLCD,
    CreateLRD,
    EnergyLQG,
    EnstrophyLQG,
    LQGParams,
    DefLQGParams,
    LQGVortex,
    DefLQGVortex,

    ## SQG

    BuildLinSysSQG,
    Calc_ψb,
    CreateModonSQG,
    Eval_ψ_SQG,
    Eval_q_SQG,
    Eval_b_SQG,
    Eval_w_SQG,
    EnergySQG,
    SQGParams,
    DefSQGParams,
    SQGVortex,
    DefSQGVortex,

    ## utils

    Calc_uv,
    Calc_∇,
    Calc_ζ

# Include all function and type definitions

## core

include("core/evp_solve.jl")
include("core/grid.jl")
include("core/num_integ.jl")
include("core/plotting.jl")
include("core/zernike.jl")

## extras

include("extras/monopoles.jl")

## LQG

include("LQG/lqg_diag.jl")
include("LQG/lqg_evp.jl")
include("LQG/lqg_high.jl")
include("LQG/lqg_low.jl")
include("LQG/lqg_type.jl")

## SQG

include("SQG/sqg_diag.jl")
include("SQG/sqg_evp.jl")
include("SQG/sqg_high.jl")
include("SQG/sqg_low.jl")
include("SQG/sqg_type.jl")

## utils

include("utils/field_utils.jl")

end

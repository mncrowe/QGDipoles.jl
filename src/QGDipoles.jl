module QGDipoles

using Jacobi
using FFTW
using LinearAlgebra
using SpecialFunctions
using QuadGK
using NLsolve

export ZernikeR, GridStruct, CreateGrid, Calc_ψq, ΔNCalc
include("create_modon.jl")

export A_func, B_func, JJ_int
include("JJ_integ.jl")

export BuildLinSys, ApplyPassiveLayers, IncludePassiveLayers, SolveInhomEVP, InhomEVP_F!, OrthogSpace 
include("lin_sys.jl")


end

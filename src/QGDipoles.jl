"""
Package for creating steady modon solutions to the layered quasi-geostrophic equations.
"""

module QGDipoles

# load required packages

using Jacobi
using FFTW
using LinearAlgebra
using SpecialFunctions
using QuadGK
using NLsolve
using CUDA

# export and define functions to calculate the numerical integrals required to build matrices

export A_func, B_func, JJ_int
include("JJ_integ.jl")

# export and define functions for building and solving the inhomogeneous eigenvalue problem

export BuildLinSys, ApplyPassiveLayers, IncludePassiveLayers, SolveInhomEVP, InhomEVP_F!, OrthogSpace 
include("lin_sys.jl")

# export and define functions to create the modon solution

export ZernikeR, GridStruct, CreateGrid, Calc_ψq, ΔNCalc, CreateModon
include("create_modon.jl")

end

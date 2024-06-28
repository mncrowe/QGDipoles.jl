
module QGDipoles

using Jacobi
using FFTW
using LinearAlgebra
using SpecialFunctions
using QuadGK

export ZernikeR
include("zernike.jl")

export A_func, B_func, JJ_int
include("JJ_integ.jl")

export BuildLinSys
include("lin_sys.jl")


end

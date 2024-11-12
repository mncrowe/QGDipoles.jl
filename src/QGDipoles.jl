"""
Package for creating steady modon solutions to the layered quasi-geostrophic equations and the
surface quasi geostrophic equations.

See `examples/` for example scripts and [here](https://mncrowe.github.io/QGDipoles.jl/) for documentation.
"""
module QGDipoles

# Load required packages

using
	Jacobi,
	FFTW,
	LinearAlgebra,
	SpecialFunctions,
	QuadGK,
	NLsolve,
	CUDA

# Export names for all functions and types

export
	# JJ_integ.jl
	A_func,
	B_func,
	JJ_int,

	# lin_sys.jl
	BuildLinSys,
	ApplyPassiveLayers,
	IncludePassiveLayers,
	SolveInhomEVP,
	InhomEVP_F!,
	OrthogSpace,

	# create_modon.jl
	ZernikeR,
	GridStruct,
	CreateGrid,
	Calc_ψq,
	Calc_ψb,
	Calc_uv,
	ΔNCalc,
	CreateModonLQG,
	CreateModonSQG,
	CreateLCD,
	CreateLRD,
	Eval_ψ_SQG,
	Eval_q_SQG,
	Eval_b_SQG,
	Eval_w_SQG,
	Calc_∇,
	CartesianGrid,
	PolarGrid

# Include all function and type definitions

include("JJ_integ.jl")
include("lin_sys.jl")
include("create_modon.jl")

end

"""
Tests for QGDipoles.jl

The tests are split across three seperate files. Each file contains testing
functions for a particular part of the code. The three files are:

- LQG_tests.jl
- SQG_tests.jl
- Wrapper_tests.jl

"""

using QGDipoles
using Test
using CUDA
using NLsolve
using SpecialFunctions

# Perform some tests on CPU if available

use_cuda = CUDA.functional() ? (false, true) : (false,)

# LQG tests

@testset "LQG_Vortices" begin

	include("LQG_tests.jl")

	# Test K matches theoretical prediction for 1-layer LCD

	@test Test1QG_K(1, 1, Inf, 0)
	
	# Test K matches semi-analytical prediction for 1-layer LRD
	
	@test Test1QG_K(1, 1, 1, 1)

	# Test K matches known numerical result for a 2-layer dipole

	@test Test2QG_K()

	for cuda_active in use_cuda	#  run for CPU & GPU if available

		# Test ψ matches known numerical result for LCD

		@test Test1QG_ψ(cuda_active)

		# Test if multi-layer PV inversion arrays are built correctly

		@test Test2QG_PVinv(cuda_active)

	end

end

# SQG tests

@testset "SQG_Vortices" begin

	include("SQG_tests.jl")

	# Test K matches known numerical result for an SQG dipole

	@test TestSQG_K()
	
	for cuda_active in use_cuda	#  run for CPU & GPU if available

		# Test maximum(v) matches known numerical solution for SQG dipole

		@test TestSQG_v(cuda_active)		# SQG vortex

	end

end

# Wrapper tests

@testset "Wrappers" begin

	include("Wrapper_tests.jl")
	
	@test true

end
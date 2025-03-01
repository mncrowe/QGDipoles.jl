"""
Tests for QGDipoles.jl

The tests are split across four seperate files. Each file contains testing
functions for a particular part of the code. The four files are:

- LQG_tests.jl
- SQG_tests.jl
- Wrapper_tests.jl
- Energetics_tests.jl
- Type_tests.jl

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

    # Test active and passive layers are correctly applied and results match
    # when layers are inculded in system or removed

    @test TestLQG_ActiveLayers()

    for cuda_active in use_cuda#  run for CPU & GPU if available

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

    for cuda_active in use_cuda#  run for CPU & GPU if available

        # Test maximum(v) matches known numerical solution for SQG dipole

        @test TestSQG_v(cuda_active)

    end

end

# Wrapper tests, use 'weird' numbers to check parameter dependence is correct

@testset "Wrappers" begin

    include("Wrapper_tests.jl")

    for cuda_active in use_cuda

        # Test wrapper for LQG, 1-layer

        @test TestWrapperLQG(1.05, 1.1, 2, 0.5; cuda = cuda_active)

        # Test wrapper for LQG, 2-layer

        @test TestWrapperLQG(0.8, 0.95, [1, 1], [0, 0.5]; cuda = cuda_active)

        # Test wrapper for SQG

        @test TestWrapperSQG(0.9, 1.1, [1, Inf], 1; cuda = cuda_active)

        # Test LCD function

        @test TestLCD(0.85, 1.05; cuda = cuda_active)

        # Test LRD function

        @test TestLRD(0.95, 0.975, 1.25, 0.45; cuda = cuda_active)

    end

end

# Energetics tests

@testset "Energetics" begin

    include("Energetics_tests.jl")

    for cuda_active in use_cuda

        # Test 1-layer LQG energy calculation

        @test TestLQGEnergy1Layer(cuda_active)

        # Test 2-layer LQG energy calculation

        @test TestLQGEnergy2Layer(cuda_active)

        # Test SQG energy calculation

        @test TestSQGEnergy(cuda_active)

    end

end

# Test examples

@testset "Examples" begin

    Ex_diag = filter(contains(r".jl$"), readdir("../examples/Diagnostics"; join=true))
    Ex_high = filter(contains(r".jl$"), readdir("../examples/High_Level"; join=true))
    Ex_lowl = filter(contains(r".jl$"), readdir("../examples/Low_Level"; join=true))
    Ex_strc = filter(contains(r".jl$"), readdir("../examples/Structures"; join=true))
    Ex = [Ex_diag; Ex_high; Ex_lowl; Ex_strc]

    for i in 1:length(Ex)

        include.(Ex[i])
        @test true

    end

end

# Test types/structures

@testset "Types" begin

    include("Type_tests.jl")

    for cuda_active in use_cuda#  run for CPU & GPU if available

        grid = CreateGrid(;
            Nx = 128,
            Ny = 128,
            Lx = [-6, 6],
            Ly = [-6, 6],
            cuda = cuda_active,
        )

        # Test LQGVortex and LQGParams construction

        @test TestLQGVortex(grid)

        # Test SQGVortex and SQGParams construction

        @test TestSQGVortex(grid)

    end

end



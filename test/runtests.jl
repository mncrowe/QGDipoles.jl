"""
Tests for QGDipoles.jl

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

    ## Test Linear System

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

        Nx, Ny = 128, 128
        Lx, Ly = 10, 10
        grid = CreateGrid(Nx, Ny, Lx, Ly; cuda = cuda_active)

        # Test ψ matches known numerical result for LCD
        @test Test1QG_ψ(grid)

        # Test if multi-layer PV inversion arrays are built correctly
        @test Test2QG_PVinv(grid)

        Nx, Ny = 256, 256
        Lx, Ly = [-10, 10], [-10, 10]
        grid = CreateGrid(Nx, Ny, Lx, Ly; cuda = cuda_active)

        # Test wrapper for LQG, 1-layer
        @test TestWrapperLQG(grid; U = 1.05, ℓ = 1.1, R = 2, β = 0.5)

        # Test wrapper for LQG, 2-layer
        @test TestWrapperLQG(grid; U = 0.8, ℓ = 0.95, R = [1, 1], β = [0, 0.5])

        # Test LCD function
        @test TestLCD(grid; U = 0.85, ℓ = 1.05)

        # Test LRD function
        @test TestLRD(grid; U = 0.95, ℓ = 0.975, R = 1.25, β = 0.45)

        # Test 1-layer LQG energy calculation
        @test TestLQGEnergy1Layer(grid)

        # Test 2-layer LQG energy calculation
        @test TestLQGEnergy2Layer(grid)

        grid = CreateGrid(;
            Nx = 128,
            Ny = 128,
            Lx = [-6, 6],
            Ly = [-6, 6],
            cuda = cuda_active,
        )

        # Test LQGVortex and LQGParams construction
        @test TestLQGVortex(grid)

    end

end

# SQG tests

@testset "SQG_Vortices" begin

    include("SQG_tests.jl")

    ## Test Linear System

    # Test K matches known numerical result for an SQG dipole
    @test TestSQG_K()

    for cuda_active in use_cuda#  run for CPU & GPU if available

        Nx, Ny = 128, 128
        Lx, Ly = 10, 10
        grid = CreateGrid(Nx, Ny, Lx, Ly; cuda = cuda_active)

        # Test maximum(v) matches known numerical solution for SQG dipole
        @test TestSQG_v(grid)

        # Test wrapper for SQG
        @test TestWrapperSQG(grid; U = 0.9, ℓ = 1.1, R = [1, Inf], β = 1)

        Nx, Ny = 256, 256
        Lx, Ly = [-10, 10], [-10, 10]
        grid = CreateGrid(Nx, Ny, Lx, Ly; cuda = cuda_active)

        # Test SQG energy calculation
        @test TestSQGEnergy(grid)

        grid = CreateGrid(;
            Nx = 128,
            Ny = 128,
            Lx = [-6, 6],
            Ly = [-6, 6],
            cuda = cuda_active,
        )

        # Test SQGVortex and SQGParams construction
        @test TestSQGVortex(grid)

    end

end

# Test examples run

@testset "Examples" begin

    Ex_dir = joinpath(pkgdir(QGDipoles), "examples")

    Ex_diag = filter(contains(r".jl$"), readdir(joinpath(Ex_dir, "Diagnostics"); join=true))
    Ex_high = filter(contains(r".jl$"), readdir(joinpath(Ex_dir, "High_Level"); join=true))
    Ex_lowl = filter(contains(r".jl$"), readdir(joinpath(Ex_dir, "Low_Level"); join=true))
    Ex_strc = filter(contains(r".jl$"), readdir(joinpath(Ex_dir, "Structures"); join=true))

    Ex = [Ex_diag; Ex_high; Ex_lowl; Ex_strc]

    for i in 1:length(Ex)

        include.(Ex[i])
        @test true

    end

end

# Test extra functions

@testset "Extras" begin

    include("extras_tests.jl")

    for cuda_active in use_cuda#  run for CPU & GPU if available

        Nx, Ny = 256, 256
        Lx, Ly = [-10, 10], [-10, 10]
        grid = CreateGrid(Nx, Ny, Lx, Ly; cuda = cuda_active)

        # Test 1 layer monopolar vortices for NaN/Inf values
        @test TestMonopoles1L(grid)

        # Test 2 layer monopolar vortices for NaN/Inf values
        @test TestMonopolesML(grid)

    end

end

nothing

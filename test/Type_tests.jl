"""
Tests for the type/structures from QGDipoles.jl.
"""

"""
    TestLQGVortex(grid)

Tests the construction of LQG vortex parameters and structures by building a vortex with
non-default parameters

"""
function TestLQGVortex(grid)

    vortex = DefLQGVortex(
        grid;
        U = 0.95,
        ℓ = 1.05,
        R = [1.1, 0.9],
        β = [0.8, 1.2],
        ActiveLayers = [1, 1],
        H = [0.25, 0.75],
        x₀ = [0.01, -0.01],
        α = 0.1,
        M = 6,
        tol = 1e-5,
        K₀ = [4.1, 4.2],
        a₀ = nothing,
        UseAnalytic = false,
        CalcVelocity = true,
        CalcVorticity = true,
        CalcEnergy = true,
        CalcEnstrophy = true,
    )

    params = UpdateParams(vortex.params; U = 1.1)

    vortex2 = UpdateVortex(grid, vortex; ℓ = 1)

    return true

end

"""
    TestSQGVortex(grid)

Tests the construction of SQG vortex parameters and structures by building a vortex with
non-default parameters

"""
function TestSQGVortex(grid)

    vortex = DefLQGVortex(
        grid;
        U = 0.95,
        ℓ = 1.05,
        R = [1.1, 0.9],
        β = 0.8,
        x₀ = [0.01, -0.01],
        α = 0.1,
        M = 10,
        tol = 1e-5,
        K₀ = 4.05,
        a₀ = nothing,
        CalcVelocity = true,
        CalcVorticity = true,
        CalcEnergy = true,
    )

    params = UpdateParams(vortex.params; U = 1.1)

    vortex2 = UpdateVortex(grid, vortex; ℓ = 1)

    return true

end


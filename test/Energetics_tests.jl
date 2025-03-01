"""
Tests for the energetics functions from QGDipoles.jl

"""

"""
    TestLQGEnergy1Layer(cuda)

Tests the values of KE and PE for a 1-layer QG model against known values.

"""
function TestLQGEnergy1Layer(cuda)

    # Set parameters

    U, ℓ, R, β = 1, 1, 1, 1
    M = 8
    tol = 1e-8

    # Create grid

    Nx, Ny = 256, 256
    Lx, Ly = 20, 20

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

    # Create LRD using numerical method

    ψ, _, _, _ = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

    KE, PE = EnergyLQG(grid, ψ; R)

    TestKE = KE[1] - 12.653032138613069 < 1e-6
    TestPE = PE[1] - 2.2180534691913154 < 1e-6

    return TestKE & TestPE

end

"""
    TestLQGEnergy2Layer(cuda)

Tests the values of KE and PE for a 2-layer QG model against known values.

"""
function TestLQGEnergy2Layer(cuda)

    # Set parameters

    U, ℓ, R, β = 1, 1, [1, 1], [0, 1]
    M = 8
    tol = 1e-8

    # Create grid

    Nx, Ny = 256, 256
    Lx, Ly = 20, 20

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

    # Create LRD using numerical method

    ψ, _, _, _ = CreateModonLQG(grid; U, ℓ, R, β, M, tol)

    KE, PE = EnergyLQG(grid, ψ; R, H = [1, 1])

    TestKE1 = KE[1] - 3.5838305725995547 < 1e-6
    TestKE2 = KE[2] - 4.7052971811344300 < 1e-6
    TestPE = PE[1] - 0.0309253993785872 < 1e-6

    return TestKE1 & TestKE2 & TestPE

end

"""
    TestSQGEnergy(cuda)

Tests the values of E and SPE for an SQG model against known values.

"""
function TestSQGEnergy(cuda)

    # Set parameters

    U, ℓ, R, β = 1, 1, [Inf, Inf], 0
    M = 8
    tol = 1e-8

    # Create grid

    Nx, Ny = 256, 256
    Lx, Ly = 20, 20

    grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

    # Create LRD using numerical method

    ψ, b, _, _ = CreateModonSQG(grid; U, ℓ, R, β, M, tol)

    E, SPE = EnergySQG(grid, ψ, b; R′ = R[2])

    TestE = E[1] - 9.750526656500448 < 1e-6
    TestSPE = SPE[1] - 30.137121991640193 < 1e-6

    return TestE & TestSPE

end

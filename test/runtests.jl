using QGDipoles
using Test

@testset "QGDipoles.jl" begin

	M, tol = 8, 1e-6

	A, B, c, d = BuildLinSys(M, 0, 0; tol)
	K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, method=1)
	
	@test abs(K[1] - 3.83171) < 1e-5

	grid = CreateGrid(128, 128, 10, 10)
	ψ, q = Calc_ψq(a, 1, 1, Inf, 0, grid)

	@test abs(maximum(ψ) - 1.28110) < 1e-5

	M, tol = 13, 1e-6

	A, B, c, d = BuildLinSys(M, [0, 0], 0; tol, sqg=true)
	K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, sqg=true)

	@test abs(K[1] - 4.12126) < 1e-5

	ψ, b = Calc_ψb(a, 1, 1, [Inf, Inf], 0, grid)
	u, v = Calc_uv(ψ, grid)

	@test abs(maximum(v) - 3.54912) < 1e-5

end

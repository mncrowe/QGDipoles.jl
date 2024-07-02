using QGDipoles
using Test

@testset "QGDipoles.jl" begin

	j₁₁ = 3.83171
	M, tol = 8, 1e-6

	A, B, c, d = BuildLinSys(M, 0, 0; tol)
	K, a = SolveInhomEVP(A, B, c, d; K₀ = [4], tol)
	
	@test abs(K[1] - 3.83171) < 1e-5

	grid = CreateGrid(128, 128, 10, 10)
	ψ, q = Calc_ψq(a, 1, 1, Inf, 0, grid)

	@test abs(maximum(ψ) - 1.2811) < 1e-5

end

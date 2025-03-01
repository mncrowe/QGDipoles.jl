"""
Test functions for functions in src/extras/

"""

"""
    TestMonopoles1L(grid)

"""
function TestMonopoles1L(grid)

    ψ, q, u, v = CreateRankine(grid; ℓ = 1.1, Γ = 1, x₀ = [0, 0])

    No_Infs_1 = ((~maximum(isinf.(ψ))) & (~maximum(isinf.(q)))) &
        ((~maximum(isinf.(u))) & (~maximum(isinf.(v))))
    No_NaNs_1 = ((~maximum(isnan.(ψ))) & (~maximum(isnan.(q)))) &
        ((~maximum(isnan.(u))) & (~maximum(isnan.(v))))

    ψ, q, u, v = Create1LMonopole(grid; ℓ = 1, Γ = 6, R = Inf, x₀ = [0.1, -0.1])

    No_Infs_2 = ((~maximum(isinf.(ψ))) & (~maximum(isinf.(q)))) &
        ((~maximum(isinf.(u))) & (~maximum(isinf.(v))))
    No_NaNs_2 = ((~maximum(isnan.(ψ))) & (~maximum(isnan.(q)))) &
        ((~maximum(isnan.(u))) & (~maximum(isnan.(v))))

    ψ, q, u, v = Create1LMonopole(grid; ℓ = 0.9, Γ = 5, R = 1, x₀ = [-0.1, 0.1])

    No_Infs_3 = ((~maximum(isinf.(ψ))) & (~maximum(isinf.(q)))) &
        ((~maximum(isinf.(u))) & (~maximum(isinf.(v))))
    No_NaNs_3 = ((~maximum(isnan.(ψ))) & (~maximum(isnan.(q)))) &
        ((~maximum(isnan.(u))) & (~maximum(isnan.(v))))

    ψ = InvertVorticity1LQG(grid, q; R = 1)

    No_Infs_4 = ~maximum(isinf.(ψ))
    No_NaNs_4 = ~maximum(isnan.(ψ))

    No_Infs = (No_Infs_1 & No_Infs_2) & (No_Infs_3 & No_Infs_4)
    No_NaNs = (No_NaNs_1 & No_NaNs_2) & (No_NaNs_3 & No_NaNs_4)

    return (No_Infs & No_NaNs)

end


"""
    TestMonopolesML(grid)

"""
function TestMonopolesML(grid)

    ψ, q = CreateLQGMonopole(grid; ℓ=1, E=[1, 1], R=[1, 1], x₀=[0, 0])

    No_Infs = (~maximum(isinf.(ψ))) & (~maximum(isinf.(q)))
    No_NaNs = (~maximum(isnan.(ψ))) & (~maximum(isnan.(q)))

    return (No_Infs & No_NaNs)

end
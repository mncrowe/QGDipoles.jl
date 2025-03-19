"""
    ZernikeR(n, x)

Define the Zernike radial function using the `jacobi` function from `SpecialFunctions`

# Arguments:
 - `n`: order, Integer
 - `x`: evaluation point, Number or Array

Note: this function is defined on ``[-1, 1]`` and is set to ``0`` for ``|x| > 1``
"""
function ZernikeR(n::Int, x::Union{Number,Array})

    y = @. (-1)^n * x * jacobi(2 * x^2 - 1, n, 0, 1) * (abs.(x) <= 1)

    return y
end

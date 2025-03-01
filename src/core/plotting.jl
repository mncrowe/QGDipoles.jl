"""
This file contains recipes for plotting solutions using `Plots.jl`

Note: Plots is not included in this package and is not needed.
"""


"""
    f(grid, F::AbstractArray; layer=1)

Recipe for plotting field `F` on `grid`

# Arguments:
 - `grid`: grid object
 - `F`: field, may have mltiple layers
 - `layer`: layer to plot (default: `1`)
"""
@recipe function f(grid, F::Array; layer = 1)

    # Set default attributes
    aspect_ratio --> 1
    colormap --> :balance
    xlims --> (minimum(grid.x), maximum(grid.x) + step(grid.x))
    ylims --> (minimum(grid.y), maximum(grid.y) + step(grid.y))
    xlabel --> "x"
    ylabel --> "y"

    # Return the plotted values
    grid.x, grid.y, transpose(@view(F[:, :, layer]))

end

"""
    f(grid, F::CuArray; layer=1)

Recipe for plotting field `F` on `grid`

# Arguments:
 - `grid`: grid object
 - `F`: field, may have mltiple layers
 - `layer`: layer to plot (default: `1`)
"""
@recipe function f(grid, F::CuArray; layer = 1)

    # Set default attributes
    aspect_ratio --> 1
    colormap --> :balance
    xlims --> (minimum(grid.x), maximum(grid.x) + step(grid.x))
    ylims --> (minimum(grid.y), maximum(grid.y) + step(grid.y))
    xlabel --> "x"
    ylabel --> "y"

    # Return the plotted values
    grid.x, grid.y, Array(transpose(@view(F[:, :, layer])))

end

"""
This file contains functions for building 2D grid structures.

These grid structures are designed to be consistent with `TwoDGrid` from
`FourierFlows.jl` hence `TwoDGrid` could also be used instead if, for
example, you wanted to use `QGDipoles.jl` to set up an initial condition
for `FourierFlows.jl` or `GeophysicalFlows.jl`.

"""

"""
    GridStruct

Stores the grid variables in physical and Fourier space

# Arguments:
 - `x`, `y`: x and y points in physical space, Ranges
 - `kr`, `l`: x and y points in Fourier space, Arrays
 - `Krsq`: `kr²+l²` in Fourier space, Array
"""
struct GridStruct
    # position ranges for x and y
    x::AbstractVector
    y::AbstractVector
    # wavenumber arrays in Fourier space
    kr::Union{Array{Float64},CuArray{Float64}}
    l::Union{Array{Float64},CuArray{Float64}}
    # K² = kr²+l² array in Fourier space
    Krsq::Union{Array{Float64},CuArray{Float64}}
end

"""
    CreateGrid(Nx, Ny, Lx, Ly; cuda=false)

Define the numerical grid as a `GridStruct`

# Arguments:
 - `Nx`, `Ny`: number of gridpoints in x and y directions, Integers
 - `Lx`, `Ly`: x and y domains, either vectors of endpoints or lengths, Vectors or Numbers
 - `cuda`: `true`; use CUDA CuArray for fields (default: `false`)
"""
function CreateGrid(
    Nx::Int,
    Ny::Int,
    Lx::Union{Number,Vector},
    Ly::Union{Number,Vector};
    cuda::Bool = false,
)

    if length(Lx) == 2

        x₀ = Lx[1]
        Lx = Lx[2] - Lx[1]

    else

        x₀ = -Lx / 2

    end

    if length(Ly) == 2

        y₀ = Ly[1]
        Ly = Ly[2] - Ly[1]

    else

        y₀ = -Ly / 2

    end

    Δx = Lx / Nx
    Δy = Ly / Ny

    x = range(x₀, step = Δx, length = Nx)
    y = range(y₀, step = Δy, length = Ny)
    kr = Array(reshape(rfftfreq(Nx, 2π / Δx), (Int(Nx / 2 + 1), 1)))
    l = Array(reshape(fftfreq(Ny, 2π / Δy), (1, Ny)))

    Krsq = @. kr^2 + l^2

    # Convert to CuArrays if using CUDA

    if cuda

        kr = CuArray(kr)
        l = CuArray(l)
        Krsq = CuArray(Krsq)

    end

    return GridStruct(x, y, kr, l, Krsq)

end

"""
    CreateGrid(; Nx=512, Ny=512, Lx=[-5,5], Ly=[-5,5], cuda=false)

Define the numerical grid as a `GridStruct` using a keyword-based method

# Arguments:
 - `Nx`, `Ny`: number of gridpoints in x and y directions, Integers
 - `Lx`, `Ly`: x and y domains, either vectors of endpoints or lengths, Vectors or Numbers
 - `cuda`: `true`; use CUDA CuArray for fields (default: `false`)
"""
CreateGrid(; Nx = 512, Ny = 512, Lx = [-5, 5], Ly = [-5, 5], cuda = false) =
    CreateGrid(Nx, Ny, Lx, Ly; cuda)

"""
    CartesianGrid(grid)

Formats the ``(x, y)`` ranges from `grid` as two-dimensional Arrays

# Arguments:
 - `grid`: grid structure containing `kr` and `l`
"""
function CartesianGrid(grid)

    x = reshape(Array(grid.x), :, 1)
    y = reshape(Array(grid.y), 1, :)

    return x, y

end

"""
    PolarGrid(x, y, x₀)

Calculates the polar coordinates from (`x`, `y`) as two-dimensional Array centred on `x₀`

# Arguments:
 - `x`, `y`: 2D Arrays for ``x`` and ``y``, created using `CartesianGrid`
 - `x₀`: Vector
"""
function PolarGrid(x, y, x₀::Vector = [0])

    r = @. sqrt((x - x₀[1])^2 + (y - x₀[2])^2)
    θ = @. atan(y - x₀[2], x - x₀[1])

    return r, θ

end

"""
Base.summary function for custom type `GridStruct`
"""
function Base.summary(g::GridStruct)

    Nx, Ny = length(g.x), length(g.y)

    if g.Krsq isa CuArray
        dev = "GPU"
    else
        dev = "CPU"
    end

    return string("Grid on ", dev, " with (Nx, Ny) = ", (Nx, Ny))
end

"""
Base.show function for custom type `GridStruct`
"""
function Base.show(io::IO, g::GridStruct)

    Nx, Ny = length(g.x), length(g.y)
    Δx, Δy = g.x[2] - g.x[1], g.y[2] - g.y[1]
    Lx, Ly = Nx * Δx, Ny * Δy

    if g.Krsq isa CuArray
        dev = "GPU"
    else
        dev = "CPU"
    end

    return print(
        io,
        "GridStruct\n",
        "  ├────────────────────── device: ",
        dev,
        "\n",
        "  ├─────────────── size (Lx, Ly): ",
        (Lx, Ly),
        "\n",
        "  ├───────── resolution (Nx, Ny): ",
        (Nx, Ny),
        "\n",
        "  ├─────── grid spacing (Δx, Δy): ",
        (Δx, Δy),
        "\n",
        "  └────────────────────── domain: x ∈ [$(g.x[1]), $(g.x[end])]",
        "\n",
        "                                  y ∈ [$(g.y[1]), $(g.y[end])]",
    )

end

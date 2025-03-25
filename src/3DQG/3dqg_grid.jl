"""
This file contains functions for building the 3D grid using a Fourier basis in the horizontal directions
and Chebyshev points (with spectral collocation) in the vertical direction.

...

"""

"""
    GridChebyshev(N, L)

Creates a grid of Chebyshev points of the second kind and a spectral collocation
differentiation matrix.

# Arguments:
 - `N`: number of gridpoints
 - `L`: vector of endpoints, ``z ∈ [L₁, L₂]``

"""

function GridChebyshev(N::Int, L::Vector)

    z = (L[1] + L[2]) / 2 .- (L[2] - L[1]) / 2 * cos.((0:(N-1))*π/(N-1))

    D = z' .- z

    D₁ = D
    D₁[D .== 0] .= 1

    A = ones(N, 1) * prod(D₁; dims = 1)

    M = A' ./ (A .* D')
    M[1:(N+1):N^2] .= sum(1 ./ D₁; dims = 1)' .- 1

    return z, M

end

"""
    GridStruct3D

Stores the 3D grid variables in physical and Fourier space

# Fields:
 - `x`, `y`: x and y points in physical space, Ranges
 - `kr`, `l`: x and y points in Fourier space, Arrays
 - `Krsq`: `kr²+l²` in Fourier space, Array
 - `z`: z points in physical space, Vector
 - `M`: spectral collocation matrix for differentiation in z, Matrix
"""
struct GridStruct3D
    # position ranges for x, y and z
    x::AbstractVector
    y::AbstractVector
    z::Vector
    # wavenumber arrays in Fourier space
    kr::Union{Array{Float64},CuArray{Float64}}
    l::Union{Array{Float64},CuArray{Float64}}
    # K² = kr²+l² array in Fourier space
    Krsq::Union{Array{Float64},CuArray{Float64}}
    # spectral collocation matrix in z
    M::Union{Array{Float64},CuArray{Float64}}
end

"""
    CreateGrid3D(Nx, Ny, Nz, Lx, Ly, Lz; cuda=false)

Define the numerical grid as a [`GridStruct3D`](@ref)

# Arguments:
 - `Nx`, `Ny`, `Nz`: number of gridpoints in x and y directions, Integers
 - `Lx`, `Ly`, `Lz`: x and y domains, either vectors of endpoints or lengths, Vectors or Numbers

# Keyword arguments:
 - `cuda`: `true`; use CUDA CuArray for fields (default: `false`)
"""
function CreateGrid3D(
    Nx::Int,
    Ny::Int,
    Nz::Int,
    Lx::Union{Number,Vector},
    Ly::Union{Number,Vector},
    Lz::Union{Number,Vector};
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

    if length(Lz) == 2

        z₀ = Lz[1]
        Lz = Lz[2] - Lz[1]

    else

        z₀ = -Lz

    end

    Δx = Lx / Nx
    Δy = Ly / Ny

    x = range(x₀, step = Δx, length = Nx)
    y = range(y₀, step = Δy, length = Ny)
    kr = Array(reshape(rfftfreq(Nx, 2π / Δx), (Int(Nx / 2 + 1), 1)))
    l = Array(reshape(fftfreq(Ny, 2π / Δy), (1, Ny)))

    Krsq = @. kr^2 + l^2

    z, M = GridChebyshev(Nz, [z₀, z₀ + Lz])

    # Convert to CuArrays if using CUDA

    if cuda

        kr = CuArray(kr)
        l = CuArray(l)
        Krsq = CuArray(Krsq)
        M = CuArray(M)

    end

    return GridStruct3D(x, y, z, kr, l, Krsq, M)

end

"""
Base.summary function for custom type [`GridStruct3D`](@ref)
"""
function Base.summary(g::GridStruct3D)

    Nx, Ny, Nz = length(g.x), length(g.y), length(g.z)

    if g.Krsq isa CuArray
        dev = "GPU"
    else
        dev = "CPU"
    end

    return string("3D Grid on ", dev, " with (Nx, Ny, Nz) = ", (Nx, Ny, Nz))
end

"""
Base.show function for custom type [`GridStruct3D`](@ref)
"""
function Base.show(io::IO, g::GridStruct3D)

    Nx, Ny, Nz = length(g.x), length(g.y), length(g.z)
    Δx, Δy = g.x[2] - g.x[1], g.y[2] - g.y[1]
    Lx, Ly, Lz = Nx * Δx, Ny * Δy, g.z[end] - g.z[1]

    if g.Krsq isa CuArray
        dev = "GPU"
    else
        dev = "CPU"
    end

    return print(
        io,
        "GridStruct3D\n",
        "  ├────────────────────── device: ",
        dev,
        "\n",
        "  ├─────────── size (Lx, Ly, Lz): ",
        (Lx, Ly, Lz),
        "\n",
        "  ├───── resolution (Nx, Ny, Nz): ",
        (Nx, Ny, Nz),
        "\n",
        "  ├─────── grid spacing (Δx, Δy): ",
        (Δx, Δy),
        "\n",
        "  └────────────────────── domain: x ∈ [$(g.x[1]), $(g.x[end])]",
        "\n",
        "                                  y ∈ [$(g.y[1]), $(g.y[end])]",
        "\n",
        "                                  z ∈ [$(g.z[1]), $(g.z[end])]",
    )

end


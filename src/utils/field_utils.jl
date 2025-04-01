"""
This file contains utilities that apply to all QG models.

Currently included are:

* `Calc_uv`: calculates 2D velocity using QG streamfunction

* `Calc_∇`: calculates the 2D gradient of a given field

* `Calc_ζ`: calculates the vertical vorticity using the QG streamfunction

If additional QG models are developed, these functions should be checked
to ensure they work with the output of these models.

"""

"""
    Calc_uv(grid, ψ)

Calculate the velocity fields from ``ψ`` using ``(u, v) = (-∂ψ/∂y, ∂ψ/∂x)``

# Arguments:
 - `grid`: grid structure containing `kr` and `l`
 - `ψ`: streamfunction, Array
"""
function Calc_uv(grid, ψ::Union{CuArray,Array})

    Nd = ndims(ψ)
    Nx, Ny = size(ψ)
    N = Int(length(ψ) / (Nx * Ny))

    # Fourier transform ψ

    ψh = rfft(reshape(ψ, Nx, Ny, N), [1, 2])

    # Calculate (u, v) in Fourier space using ∂/∂x = ik, ∂/∂y = il

    uh = -im .* grid.l .* ψh
    vh = im .* grid.kr .* ψh

    # Transform back to real space

    u = irfft(uh, Nx, [1, 2])
    v = irfft(vh, Nx, [1, 2])

    # Output arrays as 2D in SQG case for consistency

    if Nd == 2

        u = reshape(u, Nx, Ny)
        v = reshape(v, Nx, Ny)

    end

    return u, v

end

"""
    Calc_∇(grid, f)

Calculate the gradient ``∇f`` for a given field ``f``

# Arguments:
 - `grid`: grid structure containing `kr` and `l`
 - `f`: function, Array
"""
function Calc_∇(grid, f::Union{CuArray,Array})

    Nd = ndims(f)
    Nx, Ny = size(f)
    N = Int(length(f) / (Nx * Ny))

    # Fourier transform f

    fh = rfft(reshape(f, Nx, Ny, N), [1, 2])

    # Calculate (∂f/∂x, ∂f/∂y) in Fourier space using ∂/∂x = ik, ∂/∂y = il

    fxh = im .* grid.kr .* fh
    fyh = im .* grid.l .* fh

    # Transform back to real space

    fx = irfft(fxh, Nx, [1, 2])
    fy = irfft(fyh, Nx, [1, 2])

    # Output arrays as 2D in SQG case for consistency

    if Nd == 2

        fx = reshape(u, Nx, Ny)
        fy = reshape(v, Nx, Ny)

    end

    return fx, fy

end

"""
    Calc_ζ(grid, ψ)

Calculate the vertical vorticity using ``ζ = ∂v/∂x - ∂u/∂y = ∇²ψ``

# Arguments:
 - `grid`: grid structure containing `Krsq`
 - `ψ`: streamfunction, Array
"""
function Calc_ζ(grid, ψ::Union{CuArray,Array})

    Nd = ndims(ψ)
    Nx, Ny = size(ψ)
    N = Int(length(ψ) / (Nx * Ny))

    # Fourier transform ψ

    ψh = rfft(reshape(ψ, Nx, Ny, N), [1, 2])

    # Calculate ζ in Fourier space using ζ = ∇²ψ

    ζh = -grid.Krsq .* ψh

    # Transform back to real space

    ζ = irfft(ζh, Nx, [1, 2])

    # Output array as 2D in SQG case for consistency

    if Nd == 2

        ζ = reshape(ζ, Nx, Ny)

    end

    return ζ

end

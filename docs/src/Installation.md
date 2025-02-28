# Installation Instructions

Installing QGDipoles.jl is fairly straightforward and can be done using the Julia package manager.
Note that QGDipoles is not (currently) listed as a Julia package so **cannot** be installed using `] add QGDipoles`.
It is recommended to [create a new environment](https://pkgdocs.julialang.org/v1/environments/) for each project and install any new packages to that environment. 

## Installation

Installation may be done using the package manager directly by typing `]` at the Julia REPL and entering the following:
```julia
add https://github.com/mncrowe/QGDipoles.jl.git
instantiate
```
Alternatively, you may import the package manager and install by entering the following directly into the Julia REPL:
```julia
import Pkg
Pkg.add(url="https://github.com/mncrowe/QGDipoles.jl")
Pkg.instantiate()
```
This package is not compatible with versions of Julia earlier than 1.10 due to the `eachslice` function.

## Dependencies

This package requires the following dependencies:

* FFTW (v1.8.0)
* Jacobi (v0.7.0)
* LinearAlgebra
* NLsolve (v4.5.1)
* QuadGK (v2.9.4)
* SpecialFunctions (v2.4.0)
* CUDA (v5.4.3)

The specified versions are confirmed to work and earlier versions may also work.
These packages will be automatically installed with QGDipoles.jl and do not need to be added seperately.

# QGDipoles

[![Build Status](https://github.com/mncrowe/QGDipoles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mncrowe/QGDipoles.jl/actions/workflows/CI.yml?query=branch%3Amain) [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

## About

This package solves for steady, propagating dipolar solutions in a layered quasi-geostrophic (LQG) model or a surface quasi-geostrophic (SQG) model.
It is designed to work with GeophysicalFlows.jl and FourierFlows.jl.

## Installation

`QGDipoles.jl` can be installed using the Julia package manager. It is recommended to [create a new environment](https://pkgdocs.julialang.org/v1/environments/) for each project and install any new packages to that environment. 

Installation may be done using the package manager directly by typing `]` at the Julia REPL and entering the following:
```julia-repl
add https://github.com/mncrowe/QGDipoles.jl.git
instantiate
```
Alternatively, you may import the package manager and install by entering the following directly into the Julia REPL:
```julia-repl
import Pkg
Pkg.add(url="https://github.com/mncrowe/QGDipoles.jl")
Pkg.instantiate()
```
This package is not compatible with versions of Julia earlier than 1.10 due to the `eachslice` function.

## Examples and Documentation

See `examples/` for example scripts and `docs/documentation.md` for information on how to use this package and how it works. Full documetation is available [here](https://mncrowe.github.io/QGDipoles.jl).

## Files

This package contains the following files in `src/`:

* `create_modon.jl`: contains functions for calculating modon solutions on a grid using a given set of coefficients
* `lin_sys.jl`: contains functions for building and solving an inhomogeneous eigenvalue problem for the coefficients
* `JJ_int.jl`: contains numerical integration functions for calculating the terms in the inhomogeneous eigenvalue problem
* `QGDipoles.jl`: module file which loads dependencies and exports all functions

## Dependencies

This package requires the following dependencies:

* FFTW (v1.8.0)
* Jacobi (v0.7.0)
* LinearAlgebra
* NLsolve (v4.5.1)
* QuadGK (v2.9.4)
* SpecialFunctions (v2.4.0)
* CUDA (v5.4.3)

The specified versions are confirmed to work. Earlier versions may also work.

## Methodology

Dipolar vortex solutions are calculated using a method originally proposed for surface QG by [Johnson & Crowe 2023](https://doi.org/10.1017/jfm.2023.87) and extended to layered QG by [Crowe & Johnson 2024](https://doi.org/10.1017/jfm.2024.619).
An overview of this method is discussed in the [full documentation](https://mncrowe.github.io/QGDipoles.jl).
This code contains a Julia implementation of the MATLAB code included as supplementary material with [Crowe & Johnson 2024](https://doi.org/10.1017/jfm.2024.619) and also includes a version of the previously unrealeased SQG version discussed in [Johnson & Crowe 2023](https://doi.org/10.1017/jfm.2023.87) and [Crowe & Johnson 2023](https://doi.org/10.1017/jfm.2023.607).
For those interested in the original (layered QG only) implementation, it can be found [here](https://github.com/mncrowe/QGDipoles.m). 

## Contributing

Issues should be submitted [here](https://github.com/mncrowe/QGDipoles.jl/issues/) and discussions or questions [here](https://github.com/mncrowe/QGDipoles.jl/discussions).

We follow the [ColPrac guide for collaborative practices](https://github.com/SciML/ColPrac). New contributors should make sure to read that guide. Below are some additional practices we follow.

### Editing

The source files are stored in `src/`. It is recommended to open Julia with the `QGDipoles` project active by running `julia --project=.` from the root directory. Be careful not to add unnecessary packages to `Project.toml` and do not push a `Project.toml` file which contains `QGDipoles.jl` as a dependency.

We use `Documenter.jl` for creating the package documentation and `JuliaFormatter.jl` for consistent code formatting. These packages should NOT be added to the package dependencies of `QGDipoles.jl` and we recommend adding them to your base environment so they are available from any active environment.

### Formatting

This codebase is formatted using `JuliaFormatter.jl`. Formatting can be done by entering the following into the Julia REPL from the root directory:
```julia
import JuliaFormatter
JuliaFormatter.format(".")
```

### Documentation

The documentation is stored within `docs/` and can be built using `Documenter.jl` by running `docs/make.jl`. Full instructions for using `Documenter.jl` are available [here](https://documenter.juliadocs.org/stable/man/guide/). The documentation can be re-built by opening Julia using `julia --project=docs/` then entering the following into the REPL:
```julia
Pkg.develop(".")  # this adds the local copy of `QGDipoles.jl` to the `docs` environment 
include("docs/make.jl")
```
The newly built documentation will be available in `docs/build` and can be viewed by opening `docs/build/index.html` with a web browser.

### Tests

Tests are stored in `test/` and can be run using either
```julia
include("test/runtests.jl")
```
or by activating the package manager by entering `]` at the Julia REPL then entering `test`.

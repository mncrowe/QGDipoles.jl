# QGDipoles

[![Build Status](https://github.com/mncrowe/QGDipoles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mncrowe/QGDipoles.jl/actions/workflows/CI.yml?query=branch%3Amain)

## About

This package solves for steady, propagating dipolar solutions in a layered quasi-geostrophic (LQG) model or a surface quasi-geostrophic (SQG) model. It is designed to work with GeophysicalFlows.jl and FourierFlows.jl.

## Installation

To install use the Julia package manager:

```julia
julia> ]
(v1.10) pgk> add https://github.com/mncrowe/QGDipoles.jl.git
(v1.10) pgk> instantiate
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

Dipolar vortex solutions are calculated using a method originally proposed for surface QG by [Johnson & Crowe 2023](https://doi.org/10.1017/jfm.2023.87) and extended to layered QG by [Crowe & Johnson 2024](https://doi.org/10.1017/jfm.2024.619). An overview of this method is discussed in `docs/documentation.md`. This code contains a Julia implementation of the MATLAB code included as supplementary material with [Crowe & Johnson 2024](https://doi.org/10.1017/jfm.2024.619) and also includes a version of the previously unrealeased SQG version discussed in [Johnson & Crowe 2023](https://doi.org/10.1017/jfm.2023.87) and [Crowe & Johnson 2023](https://doi.org/10.1017/jfm.2023.607). For those interested in the original (layered QG only) implementation, it can be found [here](https://github.com/mncrowe/QGDipoles.m). 

 

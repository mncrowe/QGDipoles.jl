# QGDipoles

[![Build Status](https://github.com/mncrowe/QGDipoles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mncrowe/QGDipoles.jl/actions/workflows/CI.yml?query=branch%3Amain)

## About

This package solves for steady, propagating dipolar solutions in a layered quasi-geostrophic model or a surface quasi-geostrophic model. It is designed to work with GeophysicalFlows.jl and FourierFlows.jl.

## Installation

To install use the Julia package manager:

```julia
julia> ]
(v1.10) pgk> add https://github.com/mncrowe/QGDipoles.jl.git
(v1.10) pgk> instantiate
```

This package is not compatible with versions of Julia earlier than 1.10 due to the `eachslice` function.

## Examples

See `examples/` for example scripts.

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

Modon solutions are calculated using the method of [Crowe & Johnson](https://arxiv.org/abs/2404.07718).

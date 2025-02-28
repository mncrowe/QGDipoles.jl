# QGDipoles.jl

Documentation and examples for QGDipoles.jl by [Matthew N. Crowe](https://mncrowe.github.io/).

## What is QGDipoles.jl?

`QGDipoles.jl` is a Julia package which provides functions for evaluating dipolar vortex solutions in the surface quasi-geostrophic (SQG) and multi-layer quasi-geostrophic (LQG) models.
It is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics.

## What can I use QGDipoles.jl for?

The primary purpose of this package is to create dipolar vortex solutions for 


### Integration with GeophysicalFlows.jl

This package is designed to work with the `TwoDGrid` structure from [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) and [GeophysicalFlows.jl](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/)[^6].
As such, these functions may be used to define initial conditions for layered and surface quasi-geostrophic simulations which may run on either CPUs or GPUs.
However, [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) and [GeophysicalFlows.jl](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/) are NOT required to use this package as an alternative grid structure is available (created using `CreateGrid`), which uses the same field names as the eqivalent [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) function.

## How does QGDipoles.jl work?

This package is based on the semi-analytic theory of dipolar vortices derived in Johnson & Crowe 2023[^1] and Crowe & Johnson 2023[^2] for SQG solutions and Crowe & Johnson 2024[^3] for LQG solutions.
Details of the numerical method and problem parameters can be found on the Methodology page.
This code consists of an updated version of the MATLAB code released as supplementary material with Crowe & Johnson 2024[^3] and incorporates (unreleased) functions for the SQG problem.
For those interested in the original (LQG only) implementation, it can be found [here](https://github.com/mncrowe/QGDipoles.m).

[^1]: [Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2](https://doi.org/10.1017/jfm.2023.87).
[^2]: [Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10](https://doi.org/10.1017/jfm.2023.607).
[^3]: [Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1](https://doi.org/10.1017/jfm.2024.619).
# QGDipoles.jl

Documentation and examples for QGDipoles.jl by [Matthew N. Crowe](https://mncrowe.github.io/).

## What is QGDipoles.jl?

QGDipoles.jl is a Julia package which provides a series of functions that allow users to easily generate dipolar vortex solutions to two of the most commonly used QG models; the multi-layer QG (LQG) system and the surface QG (SQG) system.
Solutions for the three-dimensional quasi-geostrophic model and magneto-quasi-geostrophic model will be implemented in the future.

QGDipoles.jl is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics.
It may be used for semi-analytical studies of steady vortex solutions or to initialise a simulation with a steadily propagating, balanced vortex.
These vortex solutions may be used as initial conditions in full primitive equation models as well as QG models since they exist in QG balance so will not generate strong initial transients.

QGDipoles.jl is designed to be consistent with the framework of [GeophysicalFlows.jl](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/)[^4], a Julia package that contains modules for running LQG and SQG simulations on CPUs or GPUs.
As such, QGDipoles.jl accepts grid inputs generated using the `TwoDGrid` function from [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) and can generate solution arrays on both CPUs and GPUs using CUDA.jl.
However, [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) and [GeophysicalFlows.jl](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/) are NOT required to use this package as an alternative grid structure is available (created using `CreateGrid`), which uses the same field names as the eqivalent [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) function.

Full documentation exists with examples covering a range of LQG and SQG solutions and demonstrating how QGDipoles.jl may be integrated with [GeophysicalFlows.jl](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/).

## How does QGDipoles.jl work?

This package is based on the semi-analytic theory of dipolar vortices derived in Johnson & Crowe 2023[^1] and Crowe & Johnson 2023[^2] for SQG solutions and Crowe & Johnson 2024[^3] for LQG solutions.
Details of the numerical method and problem parameters can be found on the Methodology page.
This code consists of an updated version of the MATLAB code released as supplementary material with Crowe & Johnson 2024[^3] and incorporates (unreleased) functions for the SQG problem.
For those interested in the original (LQG only) implementation, it can be found [here](https://github.com/mncrowe/QGDipoles.m).

[^1]: [Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2](https://doi.org/10.1017/jfm.2023.87).
[^2]: [Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10](https://doi.org/10.1017/jfm.2023.607).
[^3]: [Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1](https://doi.org/10.1017/jfm.2024.619).
[^4]: [Constantinou et al., 2021, GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs, JOSS, 6(60), 3053](https://joss.theoj.org/papers/10.21105/joss.03053).
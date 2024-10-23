# QGDipoles.jl

Documentation and examples for QGDipoles.jl by Matthew N. Crowe.

## About

This Julia package provides functions for evaluating dipolar vortex solutions in the surface quasi-geostrophic (SQG) and multi-layer quasi-geostrophic (LQG) models. It is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics. This package is based on the semi-analytic theory of dipolar vortices derived in [1] and [2] for SQG solutions and [3] for LQG solutions. The method used a basis of orthogonal polynomials (Zernike radial functions) to convert a steady PDE into a linear algebra system which is solved using standard methods.

## Method Summary

The full method is outlined in [1], [2] and [3]. A summary is presented here such that the notation and examples presented later make some sense. We consider a dipolar vortex of radius $a$, moving with speed $U$. The streamfunction describing the flow is denoted by $\psi$ and potential vorticity (PV) anomaly by $q$. Velocities may be derived as $(u, v) = (-\partial_y\psi, \partial_x\psi)$. The streamfunction, $\psi$, and PV anomaly, $q$, are related through PV inversion. In the case of multiple layers, $\psi$ and $q$ are vector valued functions of length equal to the number of layers, $N$.

### Layered Quasi-Geostrophic (LQG) Solutions



### Surface Quasi-Geostrophic (SQG) Solutions

Documentation. Work in Progress.

Outline method, and give references (original analytical theory for SQG, modified semi-analytical theory for SQG, multi-layer theory for QG). Replaces original MATLAB version available as supplementary material in [cite].

Examples; make and plot examples based on example scripts.

## References

- [1]: Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2.
- [2]: Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10.
- [3]: Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1.


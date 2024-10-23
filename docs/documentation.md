# QGDipoles.jl

Documentation and examples for QGDipoles.jl by Matthew N. Crowe.

## About

This Julia package provides functions for evaluating dipolar vortex solutions in the surface quasi-geostrophic (SQG) and multi-layer quasi-geostrophic (LQG) models. It is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics. This package is based on the semi-analytic theory of dipolar vortices derived in [1] and [2] for SQG solutions and [3] for LQG solutions. The method used a basis of orthogonal polynomials (Zernike radial functions) to convert a steady PDE into a linear algebra system which is solved using standard methods. This code consists of an updated version of the MATLAB code released as supplementary material with [3] and incorporates (unreleased) functions for the SQG problem.

## Method Summary

The full method is outlined in [1], [2] and [3]. A summary is presented here such that the notation and examples presented later make some sense. We consider a dipolar vortex of radius $\ell$, moving with speed $U$. The streamfunction describing the flow is denoted by $\psi$ and potential vorticity (PV) anomaly by $q$. Velocities may be derived as $(u, v) = (-\partial_y\psi, \partial_x\psi)$. The streamfunction, $\psi$, and PV anomaly, $q$, are related through PV inversion. In the case of multiple layers, $\psi$ and $q$ are vector valued functions of length equal to the number of layers, $N$.

### Layered Quasi-Geostrophic (LQG) Solutions

In the LQG model, steady, propagating, dipolar vortices satisfy the relation

$$ q_i + \beta_i y = F_i(\psi_i + Uy),$$

where $\beta_i$ denotes the background PV gradient, $i \in [1,\dots,N]$ is the layer index and $F_i$ is an arbitrary (piecewise continuous) function. To proceed, we assume that $F_i(z) = (\beta_i/U) z$ for $x^2 + y^2 > \ell^2$ (outside the vortex) and $F_i(z) = -(K_i^2/\ell^2) z$ for $x^2 + y^2 < \ell^2$ (inside the vortex). Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system

$$ \left[ \textbf{A} -  \sum _{n = 1}^N K_n^2\textbf{B}_n \right] \textbf{a} = \textbf{c}_0 + \sum _{n = 1}^N K_n^2 \textbf{c}_n, $$

where $\textbf{A}$ and $\textbf{B}_n$ are matrices, $\textbf{a}$ is a vector containing the coefficients in the polynomial expansion, $\textbf{c}_j$ are vectors and the $K_n$ are defined in $F_i$ above and appear as unknown eigenvalues in the linear problem. In order to solve the system, $N$ additional conditions are required. These are $\textbf{d}_n \cdot a = 0$ for $n \in [1, \dots, N]$ where the $\textbf{d}_n$ are vectors. These conditions correspond to the requirement that the streamfunction and vorticity are continuous in each layer.

Solving this system determines the expansion coefficients and eigenvalues and hence allows $\psi$ and $q$ to be evaluated on any given spatial grid. In the one-layer case the problem reduces to known analytical solutions, such as the Lamb-Chaplygin dipole [4] and the Larichev-Reznik dipole [5].

### Surface Quasi-Geostrophic (SQG) Solutions

In the SQG model, steady, propagating, dipolar vortices satisfy the relation

$$ \left[\partial_z + \frac{1}{R'}\right] \psi = F(\psi + Uy),$$

where $\partial_z = [-\nabla^2 + \beta/U]^{1/2} \tanh [R [-\nabla^2 + \beta/U]^{1/2}]$ is a Dirichlet-Neumann operator linking the surface streamfunction, $\psi$, and the surface buoyancy, $b = \partial_z \psi$, and $(R, R')$ describes the baroclinic and barotropic Rossby radius respectively. 




## Examples

Examples; make and plot examples based on example scripts.

### Example 1: stuff...

### Example 2: more stuff ...

### Example 3: and so on ...

## Appendix

List of all parameters and definitions. List of functions.

## References

- [1]: Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2.
- [2]: Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10.
- [3]: Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1.
- [4]: LAMB, H. 1932 Hydrodynamics. Cambridge University Press.
- [5]: Larichev, V.D. & Reznik, G.M. 1976 Two-dimensional solitary Rossby waves. Dokl. Akad. Nauk SSSR, 12–13.


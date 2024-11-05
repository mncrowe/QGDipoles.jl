# QGDipoles.jl

Documentation and examples for QGDipoles.jl by [Matthew N. Crowe](https://mncrowe.github.io/).

## About

This Julia package provides functions for evaluating dipolar vortex solutions in the surface quasi-geostrophic (SQG) and multi-layer quasi-geostrophic (LQG) models.
It is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics.
This package is based on the semi-analytic theory of dipolar vortices derived in Johnson & Crowe 2023[^1] and Crowe & Johnson 2023[^2] for SQG solutions and Crowe & Johnson 2024[^3] for LQG solutions.
The method used a basis of orthogonal polynomials (Zernike radial functions) to convert a steady PDE into a linear algebra system which is solved using standard methods.
This code consists of an updated version of the MATLAB code released as supplementary material with Crowe & Johnson 2024[^3] and incorporates (unreleased) functions for the SQG problem.
For those interested in the original (LQG only) implementation, it can be found [here](https://github.com/mncrowe/QGDipoles.m).

## Method Summary

The full method is outlined in Johnson & Crowe 2023[^1], Crowe & Johnson 2023[^2] and 2024[^3].
A summary is presented here such that the notation and examples presented later make some sense.
We consider a dipolar vortex of radius ``\ell``, moving with speed ``U``.
This vortex consists of an isolated region of high vorticity with a closed streamline at ``x^2 + y^2 = \ell^2`` (in a frame co-moving with the vortex), hence fluid does not escape during propagation. Within the vortex core, ``x^2 + y^2 < \ell^2``, are two counter rotating regions, corresponding to a dipole.
The streamfunction describing the flow is denoted by ``\psi`` and potential vorticity (PV) anomaly by ``q``.
Velocities may be derived as ``(u, v) = (-\partial_y\psi, \partial_x\psi)``.
The streamfunction, ``\psi``, and PV anomaly, ``q``, are related through PV inversion. In the case of multiple layers, ``\psi`` and ``q`` are vector valued functions of length equal to the number of layers, ``N``.

### Layered Quasi-Geostrophic (LQG) Solutions

In the LQG model, steady, propagating, dipolar vortices satisfy the relation

```math
q_i + \beta_i y = F_i(\psi_i + Uy),
```

where ``\beta_i`` denotes the background PV gradient, ``i \in [1,\dots,N]`` is the layer index and ``F_i`` is an arbitrary (piecewise continuous) function.
To proceed, we assume that ``F_i(z) = (\beta_i/U) z`` for ``x^2 + y^2 > \ell^2`` (outside the vortex) and ``F_i(z) = -(K_i^2/\ell^2) z`` for ``x^2 + y^2 < \ell^2`` (inside the vortex).
Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system

```math
\left[ \textbf{A} -  \sum _{n = 1}^N K_n^2\textbf{B}_n \right] \textbf{a} = \textbf{c}_0 + \sum _{n = 1}^N K_n^2 \textbf{c}_n,
```

where ``\textbf{A}`` and ``\textbf{B}_n`` are matrices, ``\textbf{a}`` is a vector containing the coefficients in the polynomial expansion, ``\textbf{c}_j`` are vectors and the ``K_n`` are defined in ``F_i`` above and appear as unknown eigenvalues in the linear problem.
In order to solve the system, ``N`` additional conditions are required.
These are ``\textbf{d}_n \cdot textbf{a} = 0`` for ``n \in [1, \dots, N]`` where the ``\textbf{d}_n`` are vectors.
These conditions correspond to the requirement that the streamfunction and vorticity are continuous in each layer.
In principal, we have an infinite number of coefficients in ``\textbf{a}``.
However, since we know that these coefficients must decay with increasing index (since ``\psi`` is continuous), we can truncate the expansion after ``M`` terms.
The resulting linear system is of size ``MN \times MN``.

Solving this system determines the expansion coefficients and eigenvalues and hence allows ``\psi`` and ``q`` to be evaluated on any given spatial grid.
In the one-layer case the problem reduces to known analytical solutions, such as the Lamb-Chaplygin dipole[^4] and the Larichev-Reznik dipole[^5].

### Surface Quasi-Geostrophic (SQG) Solutions

In the SQG model, steady, propagating, dipolar vortices satisfy the relation

```math
\left[\partial_z + \frac{1}{R'}\right] \psi = F(\psi + Uy),
```

where

```math
\partial_z = \sqrt{-\nabla^2 + \beta/U} \hspace{5pt} \tanh \left[R \sqrt{-\nabla^2 + \beta/U} \right],
```

is a Dirichlet-Neumann operator linking the surface streamfunction, ``\psi``, and the surface buoyancy, ``b = \partial_z \psi``.
Here, ``(R, R')`` describes the baroclinic and barotropic Rossby radii and ``\beta`` is the background vorticity gradient.
We assume that ``F(z) = 0`` for ``x^2 + y^2 > \ell^2`` (outside the vortex) and ``F_i(z) = -(K/\ell) z`` for ``x^2 + y^2 < \ell^2``.
Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system

```math
\left[ \textbf{A} -  K\textbf{B} \right] \textbf{a} = \textbf{c}_0 + K \textbf{c}_1,
```

where ``\textbf{A}`` and ``\textbf{B}`` are matrices, ``\textbf{c}_i`` are vectors, ``\textbf{a}`` is a vector of coefficients and ``K`` is an eigenvalue related to ``F``.
An additional condition is required to solve this system for a unique set of ``K``.
This condition is taken to be continuity across the vortex boundary and corresponds to ``\textbf{d} \cdot \textbf{a} = 0`` for some vector ``\textbf{d}``.
In principal, we have an infinite number of coefficients in ``\textbf{a}``.
However, since we know that these coefficients must decay with increasing index (since ``\psi`` is continuous), we can truncate the expansion after ``M`` terms.
The resulting linear system is of size ``M \times M``.

Solving this linear system allows the surface streamfunction, ``\psi``, and surface bouyancy, ``b``, to be calculated.

### Solving the Linear System

Consider the multi-parameter, inhomogeneous eigenvalue problem

```math
\left[ \textbf{A} -  \sum _{n = 1}^N K_n^m\textbf{B}_n \right] \textbf{a} = \textbf{c}_0 + \sum _{n = 1}^N K_n^m \textbf{c}_n, \quad \textrm{s.t.} \quad \textbf{d}_n \cdot \textbf{a} = 0 \quad \textrm{for} \quad n \in [1, \dots N],
```

which describes both the SQG (``m, N = 1``) and LQG (``m = 2``) systems.
For ``N = 1``, this system may be converted into a quadratic eigenvalue problem and solved by standard techniques.
For ``N > 1``, existing techniques scale poorly with matrix size so we take an alternative approach and find ``(K, \textbf{a})`` using a root finding method, where the orthogonality conditions (``\textbf{d}_n \cdot \textbf{a} = 0``) are used to reduce the dimension of the space.
These two approaches are described in the Appendix of Crowe & Johnson 2024[^3].

### Recovering the Vortex Solution

Once the coefficients are determined, they are multiplied by the basis polynomials and summed on a specified numerical grid to give an intermediate field, ``\textbf{F}``.
The streamfunction, ``\psi``, potential vorticity anomaly, ``q`` (LQG), and surface buoyancy, ``b`` (SQG), are related to ``\textbf{F}`` by differential operators and may be calculated using discrete Fourier transforms.
Note that the streamfunction may decay slowly in the far-field for certain parameters so a sufficiently large domain is required to avoid Gibbs phenemenon near the domain edges.

### Integration with GeophysicalFlows.jl

This package is designed to work with the `TwoDGrid` structure from [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) and [GeophysicalFlows.jl](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/)[^6].
As such, these functions may be used to define initial conditions for layered and surface quasi-geostrophic simulations which may run on either CPUs or GPUs.
However, [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) and [GeophysicalFlows.jl](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/) are NOT required to use this package as an alternative grid structure is available (created using `CreateGrid`), which uses the same field names as the eqivalent [FourierFlows.jl](https://fourierflows.github.io/FourierFlowsDocumentation/) function.

## Equations and Parameters

This section contains a summary of the LQG and SQG systems of equations and tables which summarise all parameters used in functions and structures of QGDipoles.jl.

### The LQG System

The layered QG equations are

```math

(\partial_t - U \partial_x) q_i + J(\psi_i,q_i) + \beta_i \partial_x \psi_i = 0, \quad i \in \{1,2,\dots, N\},

```

where the potential vorticity in each layer, ``q_i``, is given in terms of the streamfunction in each layer, ``\psi_i``, by

```math

\begin{align*}
q_1 = &\, \nabla^2 \psi_1 + R_1^{-2} (\psi_2 - \psi_1),\\
q_i = &\, \nabla^2 \psi_i + R_i^{-2} (\psi_{i-1}-2\,\psi_i+\psi_{i+1}), \quad i \in \{2,\dots N-1\},\\
q_N = &\, \nabla^2 \psi_N  + R_N^{-2} (\psi_{N-1}-\psi_N),
\end{align*}

```

and ``R_i = \sqrt{g'H_i}/f`` is the Rossby radius of deformation in each layer.
The full list of parameters for the LQG system is given in the table below.

| Parameter | Description | Definition |
| :-- | :-- | :-- |
| ``U`` | vortex speed | - |
| ``\ell`` | vortex radius | - |
| ``\beta`` | background (y) vorticity gradient in each layer | - |
| ``R`` | Rossby radius in each layer | ``R_i = \sqrt {g'H_i} / f`` |
| ``\lambda`` | ratio of radius to ``R`` in each layer | ``\lambda_i = \ell / R_i`` |
| ``\mu`` | rescaled vorticity gradient | ``\mu_i = \beta_i \ell^2 / U`` |
| ``\alpha`` | angle of vortex propagation | - |
| ``x_0`` | position of vortex center | - |
| ``N`` | number of layers | - |
| ``M`` | number of terms in polynomial expansion | - |
| ``g'`` | buoyancy difference between each layer | - |
| ``f`` | Coriolis parameters | - |
| ``H`` | layer depth in each layer | - |

### The SQG System

The 3D QG equations are

```math

(\partial_t - U \partial_x) q + J(\psi,q) + \beta \partial_x \psi = 0, \quad \textrm{for} \quad z \in [-R, 0],

```

where

```math

q = \left[\frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2}\right] \psi, \quad \textrm{for} \quad z \in [-R, 0],

```

and ``R = NH/f`` is the Baroclinic Rossby radius.
Note that we have rescaled ``z`` by ``N/f`` so ``q`` and ``\psi`` are related by the 3D Laplacian operator.
The top boundary condition is taken to be

```math

(\partial_t - U \partial_x) [b + N^2 \eta] + J(\psi, b + N^2 \eta) = 0, \quad \textrm{on} \quad z  = 0,

```

and we assume that ``b = 0`` on the bottom surface, ``z = -R``.
Here, ``b = N \partial\psi/\partial z`` is the buoyancy and ``\eta = f\psi/g`` is the surface elevation.

The SQG system is typically derived by assuming that ``q = 0`` in the interior.
We instead take ``q = (\beta/U)\psi`` which satisfies the steady evolution equation for ``q`` given above and reduces to the usual result for ``\beta = 0``.
Since ``b(z = 0)`` can be determined from ``\psi(z = 0)`` using the Dirichlet-Neumann operator given above, this system reduces to a 2D system for the modified surface buoyancy, ``b + N^2 \eta``, only.

QGDipoles.jl solves for steady, dipolar solutions to this surface equation and hence calculates only the surface values of ``b`` and ``\psi`` using `Calc_ψb` or `CreateModonSQG`.
If a 3D solution is required, the functions `Eval_ψ_SQG`, `Eval_q_SQG` and `Eval_b_SQG` can be used to calculate `ψ`, `q` and `b` at specified depths.
Alternatively, a layered model with a large number of layers can be used to model the continuous system.
The full list of parameters for the SQG system is given in the table below.

Note: this package returns ``b/N`` rather than ``b``. When working with dimensional variables, this factor of ``1/N`` should be included manually.

| Parameter | Description | Definition |
| :-- | :-- | :-- |
| ``U`` | vortex speed | - |
| ``\ell`` | vortex radius | - |
| ``\beta`` | background (y) vorticity gradient in each layer | - |
| ``R`` | baroclinic Rossby radius | ``R = NH / f`` |
| ``R'`` | reduced barotropic Rossby radius | ``R' = R_0^2 / R`` |
| ``\lambda`` | ratio of radius to ``R`` | ``\lambda = \ell / R`` |
| ``\mu`` | rescaled vorticity gradient | ``\mu = \beta \ell^2 / U`` |
| ``\alpha`` | angle of vortex propagation | - |
| ``x_0`` | position of vortex center | - |
| ``M`` | number of terms in polynomial expansion | - |
| ``N`` | buoyancy frequency | - |
| ``R_0`` | barotropic Rossby radius | ``R_0 = \sqrt {gH} / f`` |
| ``g`` | gravitational acceleration | - |
| ``f`` | Coriolis parameters | - |
| ``H`` | layer depth | - |

## Index of Functions

This list contains all Functions and Structures from the QGDipoles.jl package.
Full descriptions can be found on the [List of Functions](https://mncrowe.github.io/QGDipoles.jl/dev/Functions/) page.

| Name | Location | Type | Description |
| :--| :-- | :-- | :-- |
| `A_func` | `src/JJ_integ.jl` | Function | Evaluates a function required to calculate the matrix ``\textbf{A}`` in the LQG case |
| `B_func` | `src/JJ_integ.jl` | Function | Evaluates a function required to calculate the matrix ``\textbf{B}`` in the LQG case |
| `JJ_int` | `src/JJ_integ.jl` | Function | Calculates a double Bessel function integral required to calculate ``\textbf{A}`` and ``\textbf{B}`` |
| `BuildLinSys` | `src/lin_sys.jl` | Function | Builds the terms in the inhomogeneous eigenvalue problem; ``\textbf{A}``, ``\textbf{B}``, ``\textbf{c}`` and ``\textbf{d}`` |
| `ApplyPassiveLayers` | `src/lin_sys.jl` | Function | Removes rows and columns corresponding to passive layers from the linear system |
|`IncludePassiveLayers` | `src/lin_sys.jl` | Function | Includes columns corresponding to passive layers in the eigenvalue and coefficient arrays |
| `SolveInhomEVP` | `src/lin_sys.jl` | Function | Solves the inhomogeneous eigenvalue problem using nonlinear root finding (``N > 1``) or standard eigenvalue methods for quadratic problems (``N = 1``) |
| `InhomEVP_F!` | `src/lin_sys.jl` | Function | Calculates the function required for `SolveInhomEVP` and it's derivatives |
| `OrthogSpace ` | `src/lin_sys.jl` | Function | Extends the input to an orthonormal basis over ``R^n`` using the Gram-Schmidt method, required for `SolveInhomEVP` |
| `ZernikeR` | `src/create_modon.jl` | Function | Define the Zernike radial function using the `jacobi` function from `SpecialFunctions.jl` |
| `GridStruct` | `src/create_modon.jl` | Structure | A structure that stores the grid variables in physical and Fourier space |
| `CreateGrid` | `src/create_modon.jl` | Function | Define the numerical grid in the form of a `GridStruct` structure |
| `Calc_ψq` | `src/create_modon.jl` | Function | Calculate ``\psi`` and ``q`` in a layered QG model using coefficients and vortex parameters |
| `Calc_ψb` | `src/create_modon.jl` | Function | Calculate ``\psi`` and ``b`` in the SQG model using coefficients and vortex parameters |
| `Calc_uv` | `src/create_modon.jl` | Function | Calculate the velocity fields from ``\psi`` using ``(u, v) = (-\partial\psi/\partial y, \partial\psi/\partial x)`` |
| `ΔNCalc` | `src/create_modon.jl` | Function | Defines a matrix used to invert for ``\psi`` and ``q`` in Fourier space |
| `CreateModonLQG` | `src/create_modon.jl` | Function | High level wrapper function for calculating ``\psi``, ``q``, ``K`` and ``\textbf{a}`` for the Layered QG model using given parameters |
| `CreateModonSQG` | `src/create_modon.jl` | Function | High level wrapper function for calculating ``\psi``, ``b``, ``K`` and ``\textbf{a}`` for the SQG model using given parameters |
| `CreateLCD` | `src/create_modon.jl` | Function | High level wrapper function for calculating ``\psi``, ``q`` and ``K`` for the Lamb-Chaplygin dipole using given parameters |
| `CreateLRD` | `src/create_modon.jl` | Function | High level wrapper function for calculating ``\psi``, ``q`` and ``K`` for the Larichev-Reznik dipole using given parameters |
| `Eval_ψ_SQG` | `src/create_modon.jl` | Function | Evaluates ψ at specified depths in the SQG model |
| `Eval_q_SQG` | `src/create_modon.jl` | Function | Evaluates q at specified depths in the SQG model |
| `Eval_b_SQG` | `src/create_modon.jl` | Function | Evaluates b at specified depths in the SQG model |


[^1]: [Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2](https://doi.org/10.1017/jfm.2023.87).
[^2]: [Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10](https://doi.org/10.1017/jfm.2023.607).
[^3]: [Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1](https://doi.org/10.1017/jfm.2024.619).
[^4]: [Lamb, H., 1932, Hydrodynamics. Cambridge University Press](https://archive.org/details/hydrodynamics00lamb).
[^5]: [Larichev, V.D. & Reznik, G.M., 1976, Two-dimensional solitary Rossby waves, Dokl. Akad. Nauk SSSR, 12–13](https://www.researchgate.net/publication/248173065_Two-dimensional_solitary_Rossby_waves).
[^6]: [Constantinou et al., 2021, GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs, JOSS, 6(60), 3053](https://joss.theoj.org/papers/10.21105/joss.03053).
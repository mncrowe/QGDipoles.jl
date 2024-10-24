# QGDipoles.jl

Documentation and examples for QGDipoles.jl by Matthew N. Crowe.

## About

This Julia package provides functions for evaluating dipolar vortex solutions in the surface quasi-geostrophic (SQG) and multi-layer quasi-geostrophic (LQG) models. It is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics. This package is based on the semi-analytic theory of dipolar vortices derived in [1] and [2] for SQG solutions and [3] for LQG solutions. The method used a basis of orthogonal polynomials (Zernike radial functions) to convert a steady PDE into a linear algebra system which is solved using standard methods. This code consists of an updated version of the MATLAB code released as supplementary material with [3] and incorporates (unreleased) functions for the SQG problem.

## Method Summary

The full method is outlined in [1], [2] and [3]. A summary is presented here such that the notation and examples presented later make some sense. We consider a dipolar vortex of radius $\ell$, moving with speed $U$. This vortex consists of an isolated region of high vorticity with a closed streamline at $x^2 + y^2 = \ell^2$, hence fluid does not escape during propagation. Within the vortex core, $x^2 + y^2 < \ell^2$, are two counter rotating regions, corresponding to a dipole. The streamfunction describing the flow is denoted by $\psi$ and potential vorticity (PV) anomaly by $q$. Velocities may be derived as $(u, v) = (-\partial_y\psi, \partial_x\psi)$. The streamfunction, $\psi$, and PV anomaly, $q$, are related through PV inversion. In the case of multiple layers, $\psi$ and $q$ are vector valued functions of length equal to the number of layers, $N$.

### Layered Quasi-Geostrophic (LQG) Solutions

In the LQG model, steady, propagating, dipolar vortices satisfy the relation

$$ q_i + \beta_i y = F_i(\psi_i + Uy),$$

where $\beta_i$ denotes the background PV gradient, $i \in [1,\dots,N]$ is the layer index and $F_i$ is an arbitrary (piecewise continuous) function. To proceed, we assume that $F_i(z) = (\beta_i/U) z$ for $x^2 + y^2 > \ell^2$ (outside the vortex) and $F_i(z) = -(K_i^2/\ell^2) z$ for $x^2 + y^2 < \ell^2$ (inside the vortex). Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system

$$ \left[ \textbf{A} -  \sum _{n = 1}^N K_n^2\textbf{B}_n \right] \textbf{a} = \textbf{c}_0 + \sum _{n = 1}^N K_n^2 \textbf{c}_n, $$

where $\textbf{A}$ and $\textbf{B}_n$ are matrices, $\textbf{a}$ is a vector containing the coefficients in the polynomial expansion, $\textbf{c}_j$ are vectors and the $K_n$ are defined in $F_i$ above and appear as unknown eigenvalues in the linear problem. In order to solve the system, $N$ additional conditions are required. These are $\textbf{d}_n \cdot a = 0$ for $n \in [1, \dots, N]$ where the $\textbf{d}_n$ are vectors. These conditions correspond to the requirement that the streamfunction and vorticity are continuous in each layer. In principal, we have an infinite number of coefficients in $\textbf{a}$. However, since we know that these coefficients must decay with increasing index (since $\psi$ in continuous), we can truncate the expansion after $M$ terms. The resulting linear system is of size $MN \times MN$.

Solving this system determines the expansion coefficients and eigenvalues and hence allows $\psi$ and $q$ to be evaluated on any given spatial grid. In the one-layer case the problem reduces to known analytical solutions, such as the Lamb-Chaplygin dipole [4] and the Larichev-Reznik dipole [5].

### Surface Quasi-Geostrophic (SQG) Solutions

In the SQG model, steady, propagating, dipolar vortices satisfy the relation

$$ \left[\partial_z + \frac{1}{R'}\right] \psi = F(\psi + Uy),$$

where

$$\partial_z = \sqrt{-\nabla^2 + \beta/U} \hspace{5pt} \tanh \left[R \sqrt{-\nabla^2 + \beta/U} \right],$$

is a Dirichlet-Neumann operator linking the surface streamfunction, $\psi$, and the surface buoyancy, $b = \partial_z \psi$, $(R, R')$ describes the baroclinic and barotropic Rossby radii and $\beta$ in the background vorticity gradient. We assume that $F(z) = 0$ for $x^2 + y^2 > \ell^2$ (outside the vortex) and $F_i(z) = -(K/\ell) z$ for $x^2 + y^2 < \ell^2$. Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system

$$ \left[ \textbf{A} -  K\textbf{B} \right] \textbf{a} = \textbf{c}_0 + K \textbf{c}_1, $$

where $\textbf{A}$ and $\textbf{B}$ are matrices, $\textbf{c}_i$ are vectors, $\textbf{a}$ is a vector of coefficients and $K$ is an eigenvalue related to $F$. An additional condition is required to solve this system for a unique set of $K$. This condition is taken to be continuity across the vortex boundary and corresponds to $\textbf{d} \cdot a = 0$ for some vector $\textbf{d}$. In principal, we have an infinite number of coefficients in $\textbf{a}$. However, since we know that these coefficients must decay with increasing index (since $\psi$ in continuous), we can truncate the expansion after $M$ terms. The resulting linear system is of size $M \times M$.

Solving this linear system allows the surface streamfunction, $\psi$, and surface bouyancy, $b$, to be calculated.

### Solving the Linear System

Consider the multi-parameter, inhomogeneous eigenvalue problem

$$ \left[ \textbf{A} -  \sum _{n = 1}^N K_n^m\textbf{B}_n \right] \textbf{a} = \textbf{c}_0 + \sum _{n = 1}^N K_n^m \textbf{c}_n, \quad \textrm{s.t.} \quad \textbf{d}_n \cdot a = 0 \quad \textrm{for} \quad n \in [1, \dots N], $$

which describes both the SQG ($m, N = 1$) and LQG ($m = 2$) systems. For $N = 1$, this system may be converted into a quadratic eigenvalue problem and solved by standard techniques. For $N > 1$, existing techniques scale poorly with matrix size so we take an alternative approach and find $(K, \textbf{a})$ using a root finding method, where the orthogonality conditions ($\textbf{d}_n \cdot a = 0$) are used to reduce the dimension of the space. These two approaches are described in the Appendix of [3].

### Recovering the Vortex Solution

Once the coefficients are determined, they are multiplied by the basis polynomials and summed on a specified numerical grid to give an intermediate field, $\textbf{F}$. The streamfunction, $\psi$, potential vorticity anomaly, $q$ (LQG), and surface buoyancy, $b$ (SQG), are related to $\textbf{F}$ by differential operators and may be calculated using discrete Fourier transforms. Note that the streamfunction may decay slowly in the far-field for certain parameters so a sufficiently large domain is required to avoid Gibbs phenemenon near the domain edges.

### Integration with GeophysicalFlows.jl

This package is designed to work with the `TwoDGrid` structure from FourierFlows.jl and GeophysicalFlows.jl. As such, these functions may be used to define initial conditions for layered and surface quasi-geostrophic simulations which may run on either CPUs or GPUs. However, FourierFlows.jl and GeophysicalFlows.jl are NOT required for this package to function as an alternative grid (created using `CreateGrid`), using the same field names as FourierFlows.jl, is available.

## Examples

Here we present some examples which demonstrate the how to use this package. Further examples are available in the `examples/` directory.

### Example 1: 1-layer QG

Let's calculate and plot the Larichev-Reznik dipole (LRD). This diople exists on the $\beta$-plane in the equivalent barotropic model so we take $\beta = R = 1$ and consider a 1-layer solution ($N = 1$). We'll also assume unit radius and velocity, $\ell = U = 1$. Let's start by loading the package and defining some parameters.

```julia

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1	# vortex speed and radius
R = 1		# Rossby radius in each layer
β = 1		# background PV gradient in each layer

M = 8		# number of coefficients in Zernike expansion
tol = 1e-8	# maximum error in solution evaluation

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

```

We've taken $M = 8$ as this is generally a sufficient number of terms to get a relative error $< 10^{-6}$ in the final result. The tolerance, `tol`, is used in calculating the terms in the linear system and a value of $10^{-8}$ corresponds to approximately the same error as our chosen $M$ value. We're also going to build a grid with $512$ points in each direction and have taken the grid size to be $10$ in each direction, which is sufficient to capture the far-field decay of the vortex. We can now build the linear system and solve for the coefficients as follows:

```julia

# Build and solve linear system for coefficients

λ = ℓ / R
μ = β * ℓ^2/U

A, B, c, d = BuildLinSys(M, λ, μ; tol)
K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol)

```

The intermediate parameters, $\lambda$ and $\mu$, describe the rescaled vortex radius and PV gradient. Finally, we can define a grid and evaluate our streamfunction, PV and velocities using:

```julia

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly)
ψ, q = Calc_ψq(a, U, ℓ, R, β, grid)
u, v = Calc_uv(ψ, grid)

```

We can plot our solution using Plots.jl:

```julia

using Plots

heatmap(grid.x, grid.y, transpose(ψ[:,:,1]);
    colormap = :balance,
    aspect_ratio=1,
    xlims=(-Lx/2, Lx/2),
    ylims = (-Ly/2, Ly/2))

```

<img width="338" alt="image" src="https://github.com/user-attachments/assets/4e17f05d-4edf-4389-968a-bc3e17c85f95">

Note that we transpose $\psi$ when plotting as $x$ corresonds to the first dimension of $\psi$.

### Example 2: multi-layer QG

This example considers a 3-layer solution and introduces the concept of active and passive layers. We define an active layer to be a layer with a closed streamline at $x^2 + y^2 = \ell^2$ whereas a passive layer has no closed streamlines. Therefore, fluid within the vortex in an active layer remains trapped in the vortex. Conversely, fluid in the passive layer is not trapped in a vortex core but can still be affected through the change in layer thickness associated with the streamfunction in neighbouring layers. Passive layers have $F_i(z) = (\beta_i/U) z$ everywhere and hence have no eigenvalue, $K_i$, to solve for. Further, the coefficients within a passive layer are zero though the solution may still be non-zero due to the coefficients in neighbouring layers. Therefore, the corresponding linear system can be simplified by removing rows and columns corresponding to passive layers and solving the reduced system for the active layers only.

We'll start by defining some parameters:

```julia

using QGDipoles

# Set problem parameters

U, ℓ = 1, 1			# vortex speed and radius
R = [1, 1, 1]			# Rossby radius in each layer
β = [0, 0, 1]			# background PV gradient in each layer
ActiveLayers = [0, 1, 0]	# 1 => layer contains vortex region
x₀ = [5, 5]			# location of vortex center

M = 8				# number of coefficients in Zernike expansion
tol = 1e-8			# maximum error in solution evaluation

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = [0, 10], [0, 10]

```

We've assumed that only the middle layer is active. Therefore our solution will describe a mid-depth propagating structure. We've also taken a background PV gradient in the lower layer only, to represent, say, a topographic slope. Finally, we've taken our vortex to be centred at $[5, 5]$ and taken $x$ and $y$ to run from $0$ to $10$.

We start by building the full linear system:

```julia

# Build and solve linear system for coefficients

λ = ℓ ./ R
μ = β * ℓ^2/U

A, B, c, d = BuildLinSys(M, λ, μ; tol)

```

Next we remove the passive layers:

```julia

A, B, c, d = ApplyPassiveLayers(A, B, c, d, ActiveLayers)

```

We can now solve the reduced system and put the passive layers, which have $(K, \textbf{a}) = (0, \textbf{0})$, back in to ensure the sizes of $K$ and $\textbf{a}$ match the number of layers:

```julia

K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol)
K, a = IncludePassiveLayers(K, a, ActiveLayers)

```

Finally, we can calculate our solution

```julia

# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly)
ψ, q = Calc_ψq(a, U, ℓ, R, β, grid, x₀)

```

### Example 3: SQG

This example covers the SQG vortex and introduces grids on a GPU.

simple, include GPU

no figure

### Example 4: Wrappers



both LQG and SQG

no figure

### Example 5: A GeophysicalFlows.jl simulation

1 layer LCD evol steady propagation.

figure of start and end state

## Appendix

The tables below summarise all parameters used in functions and structures in QGDipoles.jl. This appendix also contains an index of all functions and structures.

### LQG Parameters

| Parameter | Description | Definition |
| ----------- | ----------- | ----------- |
| $U$ | vortex speed | - |
| $\ell$ | vortex radius | - |
| $\beta$ | background (y) vorticity gradient in each layer | - |
| $R$ | Rossby radius in each layer | $R_i = \sqrt {g'H_i} / f$ |
| $\lambda$ | ratio of radius to $R$ in each layer | $\lambda_i = \ell / R_i$ |
| $\mu$ | rescaled vorticity gradient | $\mu_i = \beta_i \ell^2 / U$ |
| $\alpha$ | angle of vortex propagation | - |
| $x_0$ | position of vortex center | - |
| $N$ | number of layers | - |
| $M$ | number of terms in polynomial expansion | - |
| $g'$ | buoyancy difference between each layer | - |
| $f$ | Coriolis parameters | - |
| $H$ | layer depth in each layer | - |

### SQG Parameters

| Parameter | Description | Definition |
| ----------- | ----------- | ----------- |
| $U$ | vortex speed | - |
| $\ell$ | vortex radius | - |
| $\beta$ | background (y) vorticity gradient in each layer | - |
| $R$ | baroclinic Rossby radius | $R = NH / f$ |
| $R'$ | reduced barotropic Rossby radius | $R' = R_0^2 / R$ |
| $\lambda$ | ratio of radius to $R$ | $\lambda = \ell / R$ |
| $\mu$ | rescaled vorticity gradient | $\mu = \beta \ell^2 / U$ |
| $\alpha$ | angle of vortex propagation | - |
| $x_0$ | position of vortex center | - |
| $M$ | number of terms in polynomial expansion | - |
| $N$ | buoyancy frequency | - |
| $R_0$ | barotropic Rossby radius | $R = \sqrt {gH} / f$ |
| $g$ | gravitational acceleration | - |
| $f$ | Coriolis parameters | - |
| $H$ | layer depth | - |

### Index

| Name | Type | Description |
| ----------- | ----------- | ----------- |
| ... | ... | ... |

## References

- [1]: Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2.
- [2]: Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10.
- [3]: Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1.
- [4]: LAMB, H. 1932 Hydrodynamics. Cambridge University Press.
- [5]: Larichev, V.D. & Reznik, G.M. 1976 Two-dimensional solitary Rossby waves. Dokl. Akad. Nauk SSSR, 12–13.


# QGDipoles.jl

Documentation and examples for QGDipoles.jl by Matthew N. Crowe.

## About

This Julia package provides functions for evaluating dipolar vortex solutions in the surface quasi-geostrophic (SQG) and multi-layer quasi-geostrophic (LQG) models. It is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics. This package is based on the semi-analytic theory of dipolar vortices derived in [1] and [2] for SQG solutions and [3] for LQG solutions. The method used a basis of orthogonal polynomials (Zernike radial functions) to convert a steady PDE into a linear algebra system which is solved using standard methods. This code consists of an updated version of the MATLAB code released as supplementary material with [3] and incorporates (unreleased) functions for the SQG problem.

## Method Summary

The full method is outlined in [1], [2] and [3]. A summary is presented here such that the notation and examples presented later make some sense. We consider a dipolar vortex of radius $\ell$, moving with speed $U$. This vortex consists of an isolated region of high vorticity with a closed streamline at $x^2 + y^2 = \ell^2$ (in a frame co-moving with the vortex), hence fluid does not escape during propagation. Within the vortex core, $x^2 + y^2 < \ell^2$, are two counter rotating regions, corresponding to a dipole. The streamfunction describing the flow is denoted by $\psi$ and potential vorticity (PV) anomaly by $q$. Velocities may be derived as $(u, v) = (-\partial_y\psi, \partial_x\psi)$. The streamfunction, $\psi$, and PV anomaly, $q$, are related through PV inversion. In the case of multiple layers, $\psi$ and $q$ are vector valued functions of length equal to the number of layers, $N$.

### Layered Quasi-Geostrophic (LQG) Solutions

In the LQG model, steady, propagating, dipolar vortices satisfy the relation

$$q_i + \beta_i y = F_i(\psi_i + Uy),$$

where $\beta_i$ denotes the background PV gradient, $i \in [1,\dots,N]$ is the layer index and $F_i$ is an arbitrary (piecewise continuous) function. To proceed, we assume that $F_i(z) = (\beta_i/U) z$ for $x^2 + y^2 > \ell^2$ (outside the vortex) and $F_i(z) = -(K_i^2/\ell^2) z$ for $x^2 + y^2 < \ell^2$ (inside the vortex). Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system

$$\left[ \textbf{A} -  \sum _{n = 1}^N K_n^2\textbf{B}_n \right] \textbf{a} = \textbf{c}_0 + \sum _{n = 1}^N K_n^2 \textbf{c}_n,$$

where $\textbf{A}$ and $\textbf{B}_n$ are matrices, $\textbf{a}$ is a vector containing the coefficients in the polynomial expansion, $\textbf{c}_j$ are vectors and the $K_n$ are defined in $F_i$ above and appear as unknown eigenvalues in the linear problem. In order to solve the system, $N$ additional conditions are required. These are $\textbf{d}_n \cdot a = 0$ for $n \in [1, \dots, N]$ where the $\textbf{d}_n$ are vectors. These conditions correspond to the requirement that the streamfunction and vorticity are continuous in each layer. In principal, we have an infinite number of coefficients in $\textbf{a}$. However, since we know that these coefficients must decay with increasing index (since $\psi$ in continuous), we can truncate the expansion after $M$ terms. The resulting linear system is of size $MN \times MN$.

Solving this system determines the expansion coefficients and eigenvalues and hence allows $\psi$ and $q$ to be evaluated on any given spatial grid. In the one-layer case the problem reduces to known analytical solutions, such as the Lamb-Chaplygin dipole [4] and the Larichev-Reznik dipole [5].

### Surface Quasi-Geostrophic (SQG) Solutions

In the SQG model, steady, propagating, dipolar vortices satisfy the relation

$$\left[\partial_z + \frac{1}{R'}\right] \psi = F(\psi + Uy),$$

where

$$\partial_z = \sqrt{-\nabla^2 + \beta/U} \hspace{5pt} \tanh \left[R \sqrt{-\nabla^2 + \beta/U} \right],$$

is a Dirichlet-Neumann operator linking the surface streamfunction, $\psi$, and the surface buoyancy, $b = \partial_z \psi$. Here, $(R, R')$ describes the baroclinic and barotropic Rossby radii and $\beta$ is the background vorticity gradient. We assume that $F(z) = 0$ for $x^2 + y^2 > \ell^2$ (outside the vortex) and $F_i(z) = -(K/\ell) z$ for $x^2 + y^2 < \ell^2$. Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system

$$\left[ \textbf{A} -  K\textbf{B} \right] \textbf{a} = \textbf{c}_0 + K \textbf{c}_1,$$

where $\textbf{A}$ and $\textbf{B}$ are matrices, $\textbf{c}_i$ are vectors, $\textbf{a}$ is a vector of coefficients and $K$ is an eigenvalue related to $F$. An additional condition is required to solve this system for a unique set of $K$. This condition is taken to be continuity across the vortex boundary and corresponds to $\textbf{d} \cdot a = 0$ for some vector $\textbf{d}$. In principal, we have an infinite number of coefficients in $\textbf{a}$. However, since we know that these coefficients must decay with increasing index (since $\psi$ in continuous), we can truncate the expansion after $M$ terms. The resulting linear system is of size $M \times M$.

Solving this linear system allows the surface streamfunction, $\psi$, and surface bouyancy, $b$, to be calculated.

### Solving the Linear System

Consider the multi-parameter, inhomogeneous eigenvalue problem

$$\left[ \textbf{A} -  \sum _{n = 1}^N K_n^m\textbf{B}_n \right] \textbf{a} = \textbf{c}_0 + \sum _{n = 1}^N K_n^m \textbf{c}_n, \quad \textrm{s.t.} \quad \textbf{d}_n \cdot a = 0 \quad \textrm{for} \quad n \in [1, \dots N],$$

which describes both the SQG ($m, N = 1$) and LQG ($m = 2$) systems. For $N = 1$, this system may be converted into a quadratic eigenvalue problem and solved by standard techniques. For $N > 1$, existing techniques scale poorly with matrix size so we take an alternative approach and find $(K, \textbf{a})$ using a root finding method, where the orthogonality conditions ($\textbf{d}_n \cdot a = 0$) are used to reduce the dimension of the space. These two approaches are described in the Appendix of [3].

### Recovering the Vortex Solution

Once the coefficients are determined, they are multiplied by the basis polynomials and summed on a specified numerical grid to give an intermediate field, $\textbf{F}$. The streamfunction, $\psi$, potential vorticity anomaly, $q$ (LQG), and surface buoyancy, $b$ (SQG), are related to $\textbf{F}$ by differential operators and may be calculated using discrete Fourier transforms. Note that the streamfunction may decay slowly in the far-field for certain parameters so a sufficiently large domain is required to avoid Gibbs phenemenon near the domain edges.

### Integration with GeophysicalFlows.jl

This package is designed to work with the `TwoDGrid` structure from `FourierFlows.jl` and `GeophysicalFlows.jl` [6]. As such, these functions may be used to define initial conditions for layered and surface quasi-geostrophic simulations which may run on either CPUs or GPUs. However, `FourierFlows.jl` and `GeophysicalFlows.jl` are NOT required to use this package as an alternative grid structure (created using `CreateGrid`), which uses the same field names as `FourierFlows.jl`, is available.

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
    xlims = (-Lx/2, Lx/2),
    ylims = (-Ly/2, Ly/2))

```

![image](https://github.com/mncrowe/QGDipoles.jl/blob/gh-pages/docs/src/Ex_1.png)

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

This example covers the SQG vortex and introduces grids on a GPU. We'll start by defining some parameters. There are a few changes here compared to the LQG setup. Firstly, we'll need to set the flag `sqg` to `true` as the linear system is different in the SQG case and the functions assume LQG by default. Also, despite the SQG problem having only 1-layer, we enter `R` as a 2 element vector since we need both the (reduced) barotropic and baroclinic Rossby radii, $R$ and $R'$. We'll take these as $\infty$ using `Int` and note that these functions accept infinite Rossby radii in both the SQG and LQG cases. However, $R = 0$ is not valid since the QG assumptions break down in this limit. Note that we take $M = 20$ here and in general we'll need more coefficients for the SQG problem compared to the LQG problem as they decay slower with coefficient number. This is compensated by the fact that the SQG system is faster to calculate than the LQG system.

```julia
using QGDipoles

# Set problem parameters

U, ℓ = 1, 1     	# vortex speed and radius
R = [Inf, Inf]		# Baroclinic and Barotropic Rossby radii
β = 0		    	# background PV gradient in the interior

M = 20		    	# number of coefficients in Zernike expansion
tol = 1e-6	    	# maximum error in solution evaluation
cuda = false		# use CuArrays for grid
method = 0	    	# 0; eigensolve/nlsolve, 1; nlsolve
sqg = true	    	# functions use SQG functionality

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

```

We have introduced a couple of new variables here. Firstly, `cuda` is a flag that is passed to the grid object and when set to `true` will create the grid on an available GPU. Secondly, `method` is passed to the linear system solver, `SolveInhomEVP`, and determines if root-finding is used as the default method (`method = 1`) or if the problem is solved by eigenvalue methods for the 1-layer LQG and SQG problems (`method = 0`). In general, `method = 0` should be used, but if you have a good initial guess for $K$ and $\textbf{a}$, it may be faster to use `method = 1`.

Next we can build the linear system:

```julia
# Build and solve linear system for coefficients

λ = ℓ ./ R
μ = β * ℓ^2/U

A, B, c, d = BuildLinSys(M, λ, μ; tol, sqg)
K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, method, sqg)

```

And finally we can create our solution:

```julia
# Create grid and calculate streamfunctions and vorticities

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)
ψ, b = Calc_ψb(a, U, ℓ, R, β, grid)
u, v = Calc_uv(ψ, grid)

```

### Example 4: Wrappers

While the procedure outlined in Examples 1 to 3 gives an understanding our how this method works, it is often easier for users to be able to just call a single function to get the solution they want. Therefore, this package also includes wrappers for the SQG and LQG problems. Let's start with the LQG case and define some parameters:

```julia
# Set problem parameters

U, ℓ = 1, 1			# vortex speed and radius
R = [1, 1]			# Rossby radius in each layer
β = [0, 1]			# background PV gradient in each layer
ActiveLayers = [1, 1]		# 1 => layer contains vortex region
x₀ = [0, 0]			# position of vortex center

M = 8				# number of coefficients in Zernike expansion
tol = 1e-8			# maximum error in solution evaluation
cuda = false			# use CuArrays for grid
K₀, a₀ = [4, 4], Nothing	# guesses for K and a

# create grid

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

```

Most of these have been described in previous examples, but $K_0$ and $\textbf{a}$ are new. These are the initial guesses for $K$ and $\textbf{a}$ and are not required. They can be useful when doing a parameter sweep; since values from a case with similar parameters can be used to speed up the root-finding step for the new parameter case. In general, $K = 4$ is a good guess for most LQG and SQG vortices. In principle, there are a (countably) infinite set of solutions with increasing radial mode number. The solutions we normally think of as dipolar vortices are the first mode and higher modes are generally unstable [1].

Now we have our parameters, we can get our vortex solution with a single function call:

```julia
# create modon solution

ψ, q, K, a = CreateModonLQG(grid, M, U, ℓ, R, β, ActiveLayers, x₀; K₀, a₀, tol)

```

The SQG wrapper is similar. We start by defining our paramters:

```julia
# Set problem parameters

U, ℓ = 1, 1			# vortex speed and radius
R = [Inf, Inf]			# Baroclinic and Barotropic Rossby radii
β = 0				# background PV gradient in the interior
x₀ = [0, 0]			# position of vortex center

M = 20				# number of coefficients in Zernike expansion
tol = 1e-6			# maximum error in solution evaluation
cuda = false			# use CuArrays for grid
K₀, a₀ = 8, Nothing		# guesses for K and a

# Set grid parameters

Nx, Ny = 512, 512
Lx, Ly = 10, 10

# create grid

grid = CreateGrid(Nx, Ny, Lx, Ly; cuda)

```

Note that we've used $K_0 = 8$ this time. We'll see what happens when we create and plot our solution:

```julia
# create modon solution

ψ, b, K, a = CreateModonSQG(grid, M, U, ℓ, R, β, x₀; K₀, a₀, tol)

using Plots

heatmap(grid.x, grid.y, transpose(ψ);
    colormap = :balance,
    aspect_ratio=1,
    xlims = (-Lx/2, Lx/2),
    ylims = (-Ly/2, Ly/2))

```

If we look at $K$, we find that $K \approx 7.34205$ which is not the value we'd expect for the usual dipole solution. Instead, if we look at our plot, we see that it's a different solution with a mode 2 structure in the radial direction.

![image](https://github.com/mncrowe/QGDipoles.jl/blob/gh-pages/docs/src/Ex_4.png)

In addition to these wrapper functions, the functions `CreateLCD` and `CreateLCD` implement the Lamb-Chaplygin dipole [4] and Larichev-Reznik dipole [5] directly using the analytical solution for these cases.

### Example 5: A GeophysicalFlows.jl simulation

This package is designed to be compatible with `GeophysicalFlows.jl` [6] and provide a means of generating dipolar vortex initial conditions for layered QG and surface QG simulations. Here, we'll discuss a simple example of how to setup a 1-layer simulation in `GeophyiscalFlows.jl` using the Lamb-Chaplygin dipole as the initial condition. We'll also see that, as expected, the dipole retains it's form during the evolution and hence is a steady solution in a co-moving frame. Let's begin by defining some parameters for our vortex initial condition and our numerical simulation:

```julia
using GeophysicalFlows, QGDipoles

# Define vortex parameters

U, ℓ = 1, 1

# Set numerical simulation parameters

nx, ny = 1024, 1024
Lx, Ly = 20.48, 20.48
T = 10                        # simulation stop time
Δt = 0.01                     # timestep
Nt = Int(T/Δt)				  # number of timesteps
dev = GPU()
stepper = "FilteredRK4"

```

`GeophysicalFlows.jl` allows simulations to be run on a GPU so we've set `dev = CPU()` to use this functionality. `QGDipoles.jl` will construct vortex solutions as `CuArrays` using `CUDA` when given a grid that is stored on a GPU. We can now define our problem using the `SingleLayerQG` module from `GeophysicalFlows.jl`. This problem will contain a grid (`prob.grid`) that can be passed to functions from `QGDipoles.jl` in the same manner as grids contructed using `CreateGrid`.

```julia
# Define problem using SingleLayerQG from GeophysicalFlows.jl

prob = SingleLayerQG.Problem(dev;
		nx,
		ny,
		Lx,
		Ly,
		U = -U,			# background flow so vortex remains stationary
		dt = Δt,
		stepper)

```

Here, we've used a background flow which moves in the opposite direction to the vortex and with the same magnitude, `U`. Therefore, we're working in a frame co-moving with the vortex and we expect it to remain centred on the orogin throughout the evolution. Next, we'll use `CreateLCD` to create a Lamb-Chaplygin dipole and use this as our initial condition.

```julia
# Set initial condition

_, q₀, K = CreateLCD(prob.grid, U, ℓ)
q₀ = reshape(q₀, nx, ny)		# convert from size (nx, ny, 1) to size (nx, ny)
SingleLayerQG.set_q!(prob, q₀)

# Define Energy as a diagnostic for the simulation

diags = Diagnostic(SingleLayerQG.energy, prob; nsteps=Nt, freq=Int(Nt/100))

```

We've also defined a `Diagnostic` which will save the domain-averaged energy during the simulation. Finally, we can evolve the simulation in time:

```julia
# Evolve system forward in time

stepforward!(prob, diags, Nt)
SingleLayerQG.updatevars!(prob)

```

We can plot our initial condition and solution at $t = 10.0$ using:

```julia
using Plots

heatmap(prob.grid.x, prob.grid.y, device_array(CPU())(transpose(q₀));
		colormap = :balance,
		aspect_ratio=1,
		xlims = (-Lx/2, Lx/2),
		ylims = (-Ly/2, Ly/2))

heatmap(prob.grid.x, prob.grid.y, device_array(CPU())(transpose(prob.vars.q));
		colormap = :balance,
		aspect_ratio=1,
		xlims = (-Lx/2, Lx/2),
		ylims = (-Ly/2, Ly/2))

```

Note that we need to move our fields back to the CPU prior to plotting. The two plots are shown below and are approximately identical. Therefore, we observe that the vortex remains centred at the origin. Over long times, numerical error will result in the vortex moving at a slightly different speed to `U` and hence moving away from the origin.

![image](https://github.com/mncrowe/QGDipoles.jl/blob/gh-pages/docs/src/Ex_5a.png)![image](https://github.com/mncrowe/QGDipoles.jl/blob/gh-pages/docs/src/Ex_5b.png)

See the `GeophyiscalFlows.jl` documentation [here](https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/) for more details on how to run QG simulations.

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

| Name | Location | Type | Description |
| -----------| ----------- | ----------- | ----------- |
| `A_func` | `src/JJ_integ.jl` | Function | Evaluates a function required to calculate the matrix $\textbf{A}$ in the LQG case |
| `B_func` | `src/JJ_integ.jl` | Function | Evaluates a function required to calculate the matrix $\textbf{B}$ in the LQG case |
| `JJ_Int` | `src/JJ_integ.jl` | Function | Calculates a double Bessel function integral required to calculate $\textbf{A}$ and $\textbf{B}$ |
| `BuildLinSys` | `src/lin_sys.jl` | Function | Builds the terms in the inhomogeneous eigenvalue problem; $\textbf{A}$, $\textbf{B}$, $\textbf{c}$ and $\textbf{d}$ |
| `ApplyPassiveLayers` | `src/lin_sys.jl` | Function | Removes rows and columns corresponding to passive layers from the linear system |
|`IncludePassiveLayers` | `src/lin_sys.jl` | Function | Includes columns corresponding to passive layers in the eigenvalue and coefficient arrays |
| `SolveInhomEVP` | `src/lin_sys.jl` | Function | Solves the inhomogeneous eigenvalue problem using nonlinear root finding ($N > 1$) or standard eigenvalue methods for quadratic problems ($N = 1$) |
| `InhomEVP_F!` | `src/lin_sys.jl` | Function | Calculates the function required for `SolveInhomEVP` and it's derivatives |
| `OrthogSpace ` | `src/lin_sys.jl` | Function | Extends the input to an orthonormal basis over $R^n$ using the Gram-Schmidt method, required for `SolveInhomEVP` |
| `ZernikeR` | `src/create_modon.jl` | Function | Define the Zernike radial function using the `jacobi` function from `SpecialFunctions.jl` |
| `GridStruct` | `src/create_modon.jl` | Structure | A structure that stores the grid variables in physical and Fourier space |
| `CreateGrid` | `src/create_modon.jl` | Function | Define the numerical grid in the form of a `GridStruct` structure |
| `Calc_ψq` | `src/create_modon.jl` | Function | Calculate $\psi$ and $q$ in a layered QG model using coefficients and vortex parameters |
| `Calc_ψb` | `src/create_modon.jl` | Function | Calculate $\psi$ and $b$ in the SQG model using coefficients and vortex parameters |
| `Calc_uv` | `src/create_modon.jl` | Function | Calculate the velocity fields from $\psi$ using $(u, v) = (-\partial\psi/\partial y, \partial\psi/\partial x)$ |
| `ΔNCalc` | `src/create_modon.jl` | Function | Defines a matrix used to invert for $\psi$ and $q$ in Fourier space |
| `CreateModonLQG` | `src/create_modon.jl` | Function | High level wrapper function for calculating $\psi$, $q$, $K$ and $\textbf{a}$ for the Layered QG model using given parameters |
| `CreateModonSQG` | `src/create_modon.jl` | Function | High level wrapper function for calculating $\psi$, $b$, $K$ and $\textbf{a}$ for the SQG model using given parameters |
| `CreateLCD` | `src/create_modon.jl` | Function | High level wrapper function for calculating $\psi$, $q$ and $K$ for the Lamb-Chaplygin dipole using given parameters |
| `CreateLRD` | `src/create_modon.jl` | Function | High level wrapper function for calculating $\psi$, $q$ and $K$ for the Larichev-Reznik dipole using given parameters |

## References

- [1]: Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2.
- [2]: Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10.
- [3]: Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1.
- [4]: Lamb, H., 1932, Hydrodynamics. Cambridge University Press.
- [5]: Larichev, V.D. & Reznik, G.M., 1976, Two-dimensional solitary Rossby waves, Dokl. Akad. Nauk SSSR, 12–13.
- [6]: Constantinou et al., 2021, GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs, JOSS, 6(60), 3053.
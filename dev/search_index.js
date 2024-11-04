var documenterSearchIndex = {"docs":
[{"location":"Examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Here we present some examples which demonstrate the how to use this package. Further examples are available in the examples/ directory.","category":"page"},{"location":"Examples/#Example-1:-1-layer-QG","page":"Examples","title":"Example 1: 1-layer QG","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Let's calculate and plot the Larichev-Reznik dipole (LRD). This diople exists on the beta-plane in the equivalent barotropic model so we take beta = R = 1 and consider a 1-layer solution (N = 1). We'll also assume unit radius and velocity, ell = U = 1. Let's start by loading the package and defining some parameters.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using QGDipoles\n\n# Set problem parameters\n\nU, ℓ = 1, 1\t# vortex speed and radius\nR = 1\t\t# Rossby radius in each layer\nβ = 1\t\t# background PV gradient in each layer\n\nM = 8\t\t# number of coefficients in Zernike expansion\ntol = 1e-8\t# maximum error in solution evaluation\n\n# Set grid parameters\n\nNx, Ny = 512, 512\nLx, Ly = 10, 10\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We've taken M = 8 as this is generally a sufficient number of terms to get a relative error  10^-6 in the final result. The tolerance, tol, is used in calculating the terms in the linear system and a value of 10^-8 corresponds to approximately the same error as our chosen M value. We're also going to build a grid with 512 points in each direction and have taken the grid size to be 10 in each direction, which is sufficient to capture the far-field decay of the vortex. We can now build the linear system and solve for the coefficients as follows:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Build and solve linear system for coefficients\n\nλ = ℓ / R\nμ = β * ℓ^2/U\n\nA, B, c, d = BuildLinSys(M, λ, μ; tol)\nK, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"The intermediate parameters, lambda and mu, describe the rescaled vortex radius and PV gradient. Finally, we can define a grid and evaluate our streamfunction, PV and velocities using:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Create grid and calculate streamfunctions and vorticities\n\ngrid = CreateGrid(Nx, Ny, Lx, Ly)\nψ, q = Calc_ψq(a, U, ℓ, R, β, grid)\nu, v = Calc_uv(ψ, grid)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We can plot our solution using Plots.jl:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using Plots\n\nheatmap(grid.x, grid.y, transpose(ψ[:,:,1]);\n    colormap = :balance,\n    aspect_ratio=1,\n    xlims = (-Lx/2, Lx/2),\n    ylims = (-Ly/2, Ly/2))\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: image)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Note that we transpose psi when plotting as x corresonds to the first dimension of psi.","category":"page"},{"location":"Examples/#Example-2:-multi-layer-QG","page":"Examples","title":"Example 2: multi-layer QG","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"This example considers a 3-layer solution and introduces the concept of active and passive layers. We define an active layer to be a layer with a closed streamline at x^2 + y^2 = ell^2 whereas a passive layer has no closed streamlines. Therefore, fluid within the vortex in an active layer remains trapped in the vortex. Conversely, fluid in the passive layer is not trapped in a vortex core but can still be affected through the change in layer thickness associated with the streamfunction in neighbouring layers. Passive layers have F_i(z) = (beta_iU) z everywhere and hence have no eigenvalue, K_i, to solve for. Further, the coefficients within a passive layer are zero though the solution may still be non-zero due to the coefficients in neighbouring layers. Therefore, the corresponding linear system can be simplified by removing rows and columns corresponding to passive layers and solving the reduced system for the active layers only.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We'll start by defining some parameters:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using QGDipoles\n\n# Set problem parameters\n\nU, ℓ = 1, 1\t\t\t# vortex speed and radius\nR = [1, 1, 1]\t\t\t# Rossby radius in each layer\nβ = [0, 0, 1]\t\t\t# background PV gradient in each layer\nActiveLayers = [0, 1, 0]\t# 1 => layer contains vortex region\nx₀ = [5, 5]\t\t\t# location of vortex center\n\nM = 8\t\t\t\t# number of coefficients in Zernike expansion\ntol = 1e-8\t\t\t# maximum error in solution evaluation\n\n# Set grid parameters\n\nNx, Ny = 512, 512\nLx, Ly = [0, 10], [0, 10]\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We've assumed that only the middle layer is active. Therefore our solution will describe a mid-depth propagating structure. We've also taken a background PV gradient in the lower layer only, to represent, say, a topographic slope. Finally, we've taken our vortex to be centred at 5 5 and taken x and y to run from 0 to 10.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We start by building the full linear system:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Build and solve linear system for coefficients\n\nλ = ℓ ./ R\nμ = β * ℓ^2/U\n\nA, B, c, d = BuildLinSys(M, λ, μ; tol)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Next we remove the passive layers:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"A, B, c, d = ApplyPassiveLayers(A, B, c, d, ActiveLayers)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We can now solve the reduced system and put the passive layers, which have (K textbfa) = (0 textbf0), back in to ensure the sizes of K and textbfa match the number of layers:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"K, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol)\nK, a = IncludePassiveLayers(K, a, ActiveLayers)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Finally, we can calculate our solution","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Create grid and calculate streamfunctions and vorticities\n\ngrid = CreateGrid(Nx, Ny, Lx, Ly)\nψ, q = Calc_ψq(a, U, ℓ, R, β, grid, x₀)\n","category":"page"},{"location":"Examples/#Example-3:-SQG","page":"Examples","title":"Example 3: SQG","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"This example covers the SQG vortex and introduces grids on a GPU. We'll start by defining some parameters. There are a few changes here compared to the LQG setup. Firstly, we'll need to set the flag sqg to true as the linear system is different in the SQG case and the functions assume LQG by default. Also, despite the SQG problem having only 1-layer, we enter R as a 2 element vector since we need both the (reduced) barotropic and baroclinic Rossby radii, R and R. We'll take these as infty using Int and note that these functions accept infinite Rossby radii in both the SQG and LQG cases. However, R = 0 is not valid since the QG assumptions break down in this limit. Note that we take M = 20 here and in general we'll need more coefficients for the SQG problem compared to the LQG problem as they decay slower with coefficient number. This is compensated by the fact that the SQG system is faster to calculate than the LQG system.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using QGDipoles\n\n# Set problem parameters\n\nU, ℓ = 1, 1     \t# vortex speed and radius\nR = [Inf, Inf]\t\t# Baroclinic and Barotropic Rossby radii\nβ = 0\t\t    \t# background PV gradient in the interior\n\nM = 20\t\t    \t# number of coefficients in Zernike expansion\ntol = 1e-6\t    \t# maximum error in solution evaluation\ncuda = false\t\t# use CuArrays for grid\nmethod = 0\t    \t# 0; eigensolve/nlsolve, 1; nlsolve\nsqg = true\t    \t# functions use SQG functionality\n\n# Set grid parameters\n\nNx, Ny = 512, 512\nLx, Ly = 10, 10\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We have introduced a couple of new variables here. Firstly, cuda is a flag that is passed to the grid object and when set to true will create the grid on an available GPU. Secondly, method is passed to the linear system solver, SolveInhomEVP, and determines if root-finding is used as the default method (method = 1) or if the problem is solved by eigenvalue methods for the 1-layer LQG and SQG problems (method = 0). In general, method = 0 should be used, but if you have a good initial guess for K and textbfa, it may be faster to use method = 1.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Next we can build the linear system:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Build and solve linear system for coefficients\n\nλ = ℓ ./ R\nμ = β * ℓ^2/U\n\nA, B, c, d = BuildLinSys(M, λ, μ; tol, sqg)\nK, a = SolveInhomEVP(A, B, c, d; K₀ = 4, tol, method, sqg)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"And finally we can create our solution:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Create grid and calculate streamfunctions and vorticities\n\ngrid = CreateGrid(Nx, Ny, Lx, Ly; cuda)\nψ, b = Calc_ψb(a, U, ℓ, R, β, grid)\nu, v = Calc_uv(ψ, grid)\n","category":"page"},{"location":"Examples/#Example-4:-Wrappers","page":"Examples","title":"Example 4: Wrappers","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"While the procedure outlined in Examples 1 to 3 gives an understanding our how this method works, it is often easier for users to be able to just call a single function to get the solution they want. Therefore, this package also includes wrappers for the SQG and LQG problems. Let's start with the LQG case and define some parameters:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Set problem parameters\n\nU, ℓ = 1, 1\t\t\t# vortex speed and radius\nR = [1, 1]\t\t\t# Rossby radius in each layer\nβ = [0, 1]\t\t\t# background PV gradient in each layer\nActiveLayers = [1, 1]\t\t# 1 => layer contains vortex region\nx₀ = [0, 0]\t\t\t# position of vortex center\n\nM = 8\t\t\t\t# number of coefficients in Zernike expansion\ntol = 1e-8\t\t\t# maximum error in solution evaluation\ncuda = false\t\t\t# use CuArrays for grid\nK₀, a₀ = [4, 4], Nothing\t# guesses for K and a\n\n# create grid\n\ngrid = CreateGrid(Nx, Ny, Lx, Ly; cuda)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Most of these have been described in previous examples, but K_0 and textbfa are new. These are the initial guesses for K and textbfa and are not required. They can be useful when doing a parameter sweep; since values from a case with similar parameters can be used to speed up the root-finding step for the new parameter case. In general, K = 4 is a good guess for most LQG and SQG vortices. In principle, there are a (countably) infinite set of solutions with increasing radial mode number. The solutions we normally think of as dipolar vortices are the first mode and higher modes are generally unstable [1].","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Now we have our parameters, we can get our vortex solution with a single function call:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# create modon solution\n\nψ, q, K, a = CreateModonLQG(grid, M, U, ℓ, R, β, ActiveLayers, x₀; K₀, a₀, tol)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"The SQG wrapper is similar. We start by defining our paramters:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Set problem parameters\n\nU, ℓ = 1, 1\t\t\t# vortex speed and radius\nR = [Inf, Inf]\t\t\t# Baroclinic and Barotropic Rossby radii\nβ = 0\t\t\t\t# background PV gradient in the interior\nx₀ = [0, 0]\t\t\t# position of vortex center\n\nM = 20\t\t\t\t# number of coefficients in Zernike expansion\ntol = 1e-6\t\t\t# maximum error in solution evaluation\ncuda = false\t\t\t# use CuArrays for grid\nK₀, a₀ = 8, Nothing\t\t# guesses for K and a\n\n# Set grid parameters\n\nNx, Ny = 512, 512\nLx, Ly = 10, 10\n\n# create grid\n\ngrid = CreateGrid(Nx, Ny, Lx, Ly; cuda)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Note that we've used K_0 = 8 this time. We'll see what happens when we create and plot our solution:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# create modon solution\n\nψ, b, K, a = CreateModonSQG(grid, M, U, ℓ, R, β, x₀; K₀, a₀, tol)\n\nusing Plots\n\nheatmap(grid.x, grid.y, transpose(ψ);\n    colormap = :balance,\n    aspect_ratio=1,\n    xlims = (-Lx/2, Lx/2),\n    ylims = (-Ly/2, Ly/2))\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"If we look at K, we find that K approx 734205 which is not the value we'd expect for the usual dipole solution. Instead, if we look at our plot, we see that it's a different solution with a mode 2 structure in the radial direction.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: image)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"In addition to these wrapper functions, the functions CreateLCD and CreateLCD implement the Lamb-Chaplygin dipole [4] and Larichev-Reznik dipole [5] directly using the analytical solution for these cases.","category":"page"},{"location":"Examples/#Example-5:-A-GeophysicalFlows.jl-simulation","page":"Examples","title":"Example 5: A GeophysicalFlows.jl simulation","text":"","category":"section"},{"location":"Examples/","page":"Examples","title":"Examples","text":"This package is designed to be compatible with GeophysicalFlows.jl [6] and provide a means of generating dipolar vortex initial conditions for layered QG and surface QG simulations. Here, we'll discuss a simple example of how to setup a 1-layer simulation in GeophyiscalFlows.jl using the Lamb-Chaplygin dipole as the initial condition. We'll also see that, as expected, the dipole retains it's form during the evolution and hence is a steady solution in a co-moving frame. Let's begin by defining some parameters for our vortex initial condition and our numerical simulation:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using GeophysicalFlows, QGDipoles\n\n# Define vortex parameters\n\nU, ℓ = 1, 1\n\n# Set numerical simulation parameters\n\nnx, ny = 1024, 1024\nLx, Ly = 20.48, 20.48\nT = 10                        # simulation stop time\nΔt = 0.01                     # timestep\nNt = Int(T/Δt)\t\t\t\t  # number of timesteps\ndev = GPU()\nstepper = \"FilteredRK4\"\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"GeophysicalFlows.jl allows simulations to be run on a GPU so we've set dev = CPU() to use this functionality. QGDipoles.jl will construct vortex solutions as CuArrays using CUDA when given a grid that is stored on a GPU. We can now define our problem using the SingleLayerQG module from GeophysicalFlows.jl. This problem will contain a grid (prob.grid) that can be passed to functions from QGDipoles.jl in the same manner as grids contructed using CreateGrid.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Define problem using SingleLayerQG from GeophysicalFlows.jl\n\nprob = SingleLayerQG.Problem(dev;\n\t\tnx,\n\t\tny,\n\t\tLx,\n\t\tLy,\n\t\tU = -U,\t\t\t# background flow so vortex remains stationary\n\t\tdt = Δt,\n\t\tstepper)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Here, we've used a background flow which moves in the opposite direction to the vortex and with the same magnitude, U. Therefore, we're working in a frame co-moving with the vortex and we expect it to remain centred on the orogin throughout the evolution. Next, we'll use CreateLCD to create a Lamb-Chaplygin dipole and use this as our initial condition.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Set initial condition\n\n_, q₀, K = CreateLCD(prob.grid, U, ℓ)\nq₀ = reshape(q₀, nx, ny)\t\t# convert from size (nx, ny, 1) to size (nx, ny)\nSingleLayerQG.set_q!(prob, q₀)\n\n# Define Energy as a diagnostic for the simulation\n\ndiags = Diagnostic(SingleLayerQG.energy, prob; nsteps=Nt, freq=Int(Nt/100))\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We've also defined a Diagnostic which will save the domain-averaged energy during the simulation. Finally, we can evolve the simulation in time:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"# Evolve system forward in time\n\nstepforward!(prob, diags, Nt)\nSingleLayerQG.updatevars!(prob)\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"We can plot our initial condition and solution at t = 100 using:","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"using Plots\n\nheatmap(prob.grid.x, prob.grid.y, device_array(CPU())(transpose(q₀));\n\t\tcolormap = :balance,\n\t\taspect_ratio=1,\n\t\txlims = (-Lx/2, Lx/2),\n\t\tylims = (-Ly/2, Ly/2))\n\nheatmap(prob.grid.x, prob.grid.y, device_array(CPU())(transpose(prob.vars.q));\n\t\tcolormap = :balance,\n\t\taspect_ratio=1,\n\t\txlims = (-Lx/2, Lx/2),\n\t\tylims = (-Ly/2, Ly/2))\n","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"Note that we need to move our fields back to the CPU prior to plotting. The two plots are shown below and are approximately identical. Therefore, we observe that the vortex remains centred at the origin. Over long times, numerical error will result in the vortex moving at a slightly different speed to U and hence moving away from the origin.","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"(Image: image)(Image: image)","category":"page"},{"location":"Examples/","page":"Examples","title":"Examples","text":"See the GeophyiscalFlows.jl documentation here for more details on how to run QG simulations.","category":"page"},{"location":"Functions/#Functions","page":"List of Functions","title":"Functions","text":"","category":"section"},{"location":"Functions/","page":"List of Functions","title":"List of Functions","text":"ApplyPassiveLayers","category":"page"},{"location":"Functions/#QGDipoles.ApplyPassiveLayers","page":"List of Functions","title":"QGDipoles.ApplyPassiveLayers","text":"Function: ApplyPassiveLayers\n\nRemoves rows and columns corresponding to passive layers from the system\n\nArguments:\n\nA, B, c, d: inhomogeneous eigenvalue problem terms, Arrays\nActiveLayers: vector of 1s or 0s where 1 denotes an active layer, Number or Vector\n\n\n\n\n\n","category":"function"},{"location":"#QGDipoles.jl","page":"Home","title":"QGDipoles.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation and examples for QGDipoles.jl by Matthew N. Crowe.","category":"page"},{"location":"#About","page":"Home","title":"About","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This Julia package provides functions for evaluating dipolar vortex solutions in the surface quasi-geostrophic (SQG) and multi-layer quasi-geostrophic (LQG) models. It is intended for use by those researching vortex dynamics in strongly rotating flows, in particular for researchers in physical oceanography and atmospheric dynamics. This package is based on the semi-analytic theory of dipolar vortices derived in [1] and [2] for SQG solutions and [3] for LQG solutions. The method used a basis of orthogonal polynomials (Zernike radial functions) to convert a steady PDE into a linear algebra system which is solved using standard methods. This code consists of an updated version of the MATLAB code released as supplementary material with [3] and incorporates (unreleased) functions for the SQG problem.","category":"page"},{"location":"#Method-Summary","page":"Home","title":"Method Summary","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The full method is outlined in [1], [2] and [3]. A summary is presented here such that the notation and examples presented later make some sense. We consider a dipolar vortex of radius ell, moving with speed U. This vortex consists of an isolated region of high vorticity with a closed streamline at x^2 + y^2 = ell^2 (in a frame co-moving with the vortex), hence fluid does not escape during propagation. Within the vortex core, x^2 + y^2  ell^2, are two counter rotating regions, corresponding to a dipole. The streamfunction describing the flow is denoted by psi and potential vorticity (PV) anomaly by q. Velocities may be derived as (u v) = (-partial_ypsi partial_xpsi). The streamfunction, psi, and PV anomaly, q, are related through PV inversion. In the case of multiple layers, psi and q are vector valued functions of length equal to the number of layers, N.","category":"page"},{"location":"#Layered-Quasi-Geostrophic-(LQG)-Solutions","page":"Home","title":"Layered Quasi-Geostrophic (LQG) Solutions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In the LQG model, steady, propagating, dipolar vortices satisfy the relation","category":"page"},{"location":"","page":"Home","title":"Home","text":"q_i + beta_i y = F_i(psi_i + Uy)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where beta_i denotes the background PV gradient, i in 1dotsN is the layer index and F_i is an arbitrary (piecewise continuous) function. To proceed, we assume that F_i(z) = (beta_iU) z for x^2 + y^2  ell^2 (outside the vortex) and F_i(z) = -(K_i^2ell^2) z for x^2 + y^2  ell^2 (inside the vortex). Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system","category":"page"},{"location":"","page":"Home","title":"Home","text":"left textbfA -  sum _n = 1^N K_n^2textbfB_n right textbfa = textbfc_0 + sum _n = 1^N K_n^2 textbfc_n","category":"page"},{"location":"","page":"Home","title":"Home","text":"where textbfA and textbfB_n are matrices, textbfa is a vector containing the coefficients in the polynomial expansion, textbfc_j are vectors and the K_n are defined in F_i above and appear as unknown eigenvalues in the linear problem. In order to solve the system, N additional conditions are required. These are textbfd_n cdot a = 0 for n in 1 dots N where the textbfd_n are vectors. These conditions correspond to the requirement that the streamfunction and vorticity are continuous in each layer. In principal, we have an infinite number of coefficients in textbfa. However, since we know that these coefficients must decay with increasing index (since psi in continuous), we can truncate the expansion after M terms. The resulting linear system is of size MN times MN.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Solving this system determines the expansion coefficients and eigenvalues and hence allows psi and q to be evaluated on any given spatial grid. In the one-layer case the problem reduces to known analytical solutions, such as the Lamb-Chaplygin dipole [4] and the Larichev-Reznik dipole [5].","category":"page"},{"location":"#Surface-Quasi-Geostrophic-(SQG)-Solutions","page":"Home","title":"Surface Quasi-Geostrophic (SQG) Solutions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In the SQG model, steady, propagating, dipolar vortices satisfy the relation","category":"page"},{"location":"","page":"Home","title":"Home","text":"leftpartial_z + frac1Rright psi = F(psi + Uy)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where","category":"page"},{"location":"","page":"Home","title":"Home","text":"partial_z = sqrt-nabla^2 + betaU hspace5pt tanh leftR sqrt-nabla^2 + betaU right","category":"page"},{"location":"","page":"Home","title":"Home","text":"is a Dirichlet-Neumann operator linking the surface streamfunction, psi, and the surface buoyancy, b = partial_z psi. Here, (R R) describes the baroclinic and barotropic Rossby radii and beta is the background vorticity gradient. We assume that F(z) = 0 for x^2 + y^2  ell^2 (outside the vortex) and F_i(z) = -(Kell) z for x^2 + y^2  ell^2. Using a Hankel transform and expansion in term of Zernike radial functions, the problem may be reduced to the linear algebra system","category":"page"},{"location":"","page":"Home","title":"Home","text":"left textbfA -  KtextbfB right textbfa = textbfc_0 + K textbfc_1","category":"page"},{"location":"","page":"Home","title":"Home","text":"where textbfA and textbfB are matrices, textbfc_i are vectors, textbfa is a vector of coefficients and K is an eigenvalue related to F. An additional condition is required to solve this system for a unique set of K. This condition is taken to be continuity across the vortex boundary and corresponds to textbfd cdot a = 0 for some vector textbfd. In principal, we have an infinite number of coefficients in textbfa. However, since we know that these coefficients must decay with increasing index (since psi in continuous), we can truncate the expansion after M terms. The resulting linear system is of size M times M.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Solving this linear system allows the surface streamfunction, psi, and surface bouyancy, b, to be calculated.","category":"page"},{"location":"#Solving-the-Linear-System","page":"Home","title":"Solving the Linear System","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Consider the multi-parameter, inhomogeneous eigenvalue problem","category":"page"},{"location":"","page":"Home","title":"Home","text":"left textbfA -  sum _n = 1^N K_n^mtextbfB_n right textbfa = textbfc_0 + sum _n = 1^N K_n^m textbfc_n quad textrmst quad textbfd_n cdot a = 0 quad textrmfor quad n in 1 dots N","category":"page"},{"location":"","page":"Home","title":"Home","text":"which describes both the SQG (m N = 1) and LQG (m = 2) systems. For N = 1, this system may be converted into a quadratic eigenvalue problem and solved by standard techniques. For N  1, existing techniques scale poorly with matrix size so we take an alternative approach and find (K textbfa) using a root finding method, where the orthogonality conditions (textbfd_n cdot a = 0) are used to reduce the dimension of the space. These two approaches are described in the Appendix of [3].","category":"page"},{"location":"#Recovering-the-Vortex-Solution","page":"Home","title":"Recovering the Vortex Solution","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Once the coefficients are determined, they are multiplied by the basis polynomials and summed on a specified numerical grid to give an intermediate field, textbfF. The streamfunction, psi, potential vorticity anomaly, q (LQG), and surface buoyancy, b (SQG), are related to textbfF by differential operators and may be calculated using discrete Fourier transforms. Note that the streamfunction may decay slowly in the far-field for certain parameters so a sufficiently large domain is required to avoid Gibbs phenemenon near the domain edges.","category":"page"},{"location":"#Integration-with-GeophysicalFlows.jl","page":"Home","title":"Integration with GeophysicalFlows.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is designed to work with the TwoDGrid structure from FourierFlows.jl and GeophysicalFlows.jl [6]. As such, these functions may be used to define initial conditions for layered and surface quasi-geostrophic simulations which may run on either CPUs or GPUs. However, FourierFlows.jl and GeophysicalFlows.jl are NOT required to use this package as an alternative grid structure (created using CreateGrid), which uses the same field names as FourierFlows.jl, is available.","category":"page"},{"location":"#Appendix","page":"Home","title":"Appendix","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The tables below summarise all parameters used in functions and structures in QGDipoles.jl. This appendix also contains an index of all functions and structures.","category":"page"},{"location":"#LQG-Parameters","page":"Home","title":"LQG Parameters","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Parameter Description Definition\nU vortex speed -\nell vortex radius -\nbeta background (y) vorticity gradient in each layer -\nR Rossby radius in each layer R_i = sqrt gH_i  f\nlambda ratio of radius to R in each layer lambda_i = ell  R_i\nmu rescaled vorticity gradient mu_i = beta_i ell^2  U\nalpha angle of vortex propagation -\nx_0 position of vortex center -\nN number of layers -\nM number of terms in polynomial expansion -\ng buoyancy difference between each layer -\nf Coriolis parameters -\nH layer depth in each layer -","category":"page"},{"location":"#SQG-Parameters","page":"Home","title":"SQG Parameters","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Parameter Description Definition\nU vortex speed -\nell vortex radius -\nbeta background (y) vorticity gradient in each layer -\nR baroclinic Rossby radius R = NH  f\nR reduced barotropic Rossby radius R = R_0^2  R\nlambda ratio of radius to R lambda = ell  R\nmu rescaled vorticity gradient mu = beta ell^2  U\nalpha angle of vortex propagation -\nx_0 position of vortex center -\nM number of terms in polynomial expansion -\nN buoyancy frequency -\nR_0 barotropic Rossby radius R = sqrt gH  f\ng gravitational acceleration -\nf Coriolis parameters -\nH layer depth -","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Name Location Type Description\nA_func src/JJ_integ.jl Function Evaluates a function required to calculate the matrix textbfA in the LQG case\nB_func src/JJ_integ.jl Function Evaluates a function required to calculate the matrix textbfB in the LQG case\nJJ_Int src/JJ_integ.jl Function Calculates a double Bessel function integral required to calculate textbfA and textbfB\nBuildLinSys src/lin_sys.jl Function Builds the terms in the inhomogeneous eigenvalue problem; textbfA, textbfB, textbfc and textbfd\nApplyPassiveLayers src/lin_sys.jl Function Removes rows and columns corresponding to passive layers from the linear system\nIncludePassiveLayers src/lin_sys.jl Function Includes columns corresponding to passive layers in the eigenvalue and coefficient arrays\nSolveInhomEVP src/lin_sys.jl Function Solves the inhomogeneous eigenvalue problem using nonlinear root finding (N  1) or standard eigenvalue methods for quadratic problems (N = 1)\nInhomEVP_F! src/lin_sys.jl Function Calculates the function required for SolveInhomEVP and it's derivatives\nOrthogSpace src/lin_sys.jl Function Extends the input to an orthonormal basis over R^n using the Gram-Schmidt method, required for SolveInhomEVP\nZernikeR src/create_modon.jl Function Define the Zernike radial function using the jacobi function from SpecialFunctions.jl\nGridStruct src/create_modon.jl Structure A structure that stores the grid variables in physical and Fourier space\nCreateGrid src/create_modon.jl Function Define the numerical grid in the form of a GridStruct structure\nCalc_ψq src/create_modon.jl Function Calculate psi and q in a layered QG model using coefficients and vortex parameters\nCalc_ψb src/create_modon.jl Function Calculate psi and b in the SQG model using coefficients and vortex parameters\nCalc_uv src/create_modon.jl Function Calculate the velocity fields from psi using (u v) = (-partialpsipartial y partialpsipartial x)\nΔNCalc src/create_modon.jl Function Defines a matrix used to invert for psi and q in Fourier space\nCreateModonLQG src/create_modon.jl Function High level wrapper function for calculating psi, q, K and textbfa for the Layered QG model using given parameters\nCreateModonSQG src/create_modon.jl Function High level wrapper function for calculating psi, b, K and textbfa for the SQG model using given parameters\nCreateLCD src/create_modon.jl Function High level wrapper function for calculating psi, q and K for the Lamb-Chaplygin dipole using given parameters\nCreateLRD src/create_modon.jl Function High level wrapper function for calculating psi, q and K for the Larichev-Reznik dipole using given parameters","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"[1]: Johnson, E. R., and M. N. Crowe, 2023, Oceanic dipoles in a surface quasigeostrophic model, J. Fluid Mech., 958, R2.\n[2]: Crowe, M. N., and E. R. Johnson, 2023, The evolution of surface quasi-geostrophic modons on sloping topography, J. Fluid. Mech., 970, A10.\n[3]: Crowe, M. N., and E. R. Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic model, J. Fluid. Mech., 994, R1.\n[4]: Lamb, H., 1932, Hydrodynamics. Cambridge University Press.\n[5]: Larichev, V.D. & Reznik, G.M., 1976, Two-dimensional solitary Rossby waves, Dokl. Akad. Nauk SSSR, 12–13.\n[6]: Constantinou et al., 2021, GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs, JOSS, 6(60), 3053.","category":"page"}]
}
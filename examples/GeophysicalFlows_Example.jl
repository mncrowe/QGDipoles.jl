# This example is covered in example 5 of docs/documentation.md and combines GDipoles.jl with
# GeophysicalFlows.jl to simulate a steadily propagating dipole

using GeophysicalFlows, QGDipoles

# Define vortex parameters

U, ℓ = 1, 1

# Set numerical simulation parameters

nx, ny = 1024, 1024
Lx, Ly = 20.48, 20.48
T = 10
Δt = 0.01
Nt = Int(T/Δt)				# number of timesteps
dev = GPU()
stepper = "FilteredRK4"

# Define problem using SingleLayerQG from GeophysicalFlows.jl

prob = SingleLayerQG.Problem(dev;
		nx,
		ny,
		Lx,
		Ly,
		U = -U,			# background flow so vortex remains stationary
		dt = Δt,
		stepper)

# Set initial condition

_, q₀, K = CreateLCD(prob.grid, U, ℓ)
q₀ = reshape(q₀, nx, ny)		# convert from size (nx, ny, 1) to size (nx, ny)
SingleLayerQG.set_q!(prob, q₀)

# Define Energy as a diagnostic for the simulation

diags = Diagnostic(SingleLayerQG.energy, prob; nsteps=Nt, freq=Int(Nt/100))

# Evolve system forward in time

stepforward!(prob, diags, Nt)
SingleLayerQG.updatevars!(prob)

# Plot initial and final fields

using Plots

heatmap(prob.grid.x, prob.grid.y, device_array(CPU())(transpose(q₀));
		colormap = :balance,
		aspect_ratio=1,
		xlims= (-Lx/2, Lx/2),
		ylims = (-Ly/2, Ly/2))

heatmap(prob.grid.x, prob.grid.y, device_array(CPU())(transpose(prob.vars.q));
		colormap = :balance,
		aspect_ratio=1,
		xlims= (-Lx/2, Lx/2),
		ylims = (-Ly/2, Ly/2))
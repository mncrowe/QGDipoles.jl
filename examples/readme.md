# Examples

This directory contains a series of examples for generating dipolar vortex solutions in the layered QG and surface QG model. Further discussion of these examples is given [here](https://mncrowe.github.io/QGDipoles.jl/dev/Examples/).

These examples are divided across 5 directories:

* Diagnostics: examples demonstrating how to use diagnostics such as energy and enstrophy

* GeophysicalFlows: examples showing how to combine QGDipoles.jl with GeophysicalFlows.jl

* High_Level: examples showing high-level functions for calculating dipolar vortices

* Low_Level: examples showing low-level functions for calculating dipolar vortices

* Structures: examples showing functions for generating vortex solutions as structures

Note: using the low-level functions will reqire an understanding of the [methodology](https://mncrowe.github.io/QGDipoles.jl/dev/Methodology/).

# Available Scripts

The available example scripts are:

### Diagnostics

* Energy_LQG.jl: Constructs a 2-layer dipolar vortex solution and calculates the kinetic energy, potential energy and enstrophy

* Energy_SQG.jl: Constructs an SQG dipolar vortex solution and calculates the domain integrated energy and the surface potential energy

### GeophysicalFlows

* GeophysicalFlows_Example.jl: Demostrates how to used `QGDipoles.jl` to generate initial conditions for `GeophysicalFlows.jl` (corresponds to Example 5 [here](https://mncrowe.github.io/QGDipoles.jl/dev/Examples/)).

### High_Level

* Wrapper_LQG.jl: Demostrates the use of wrapper functions to easily create layered QG solutions (corresponds to Example 4 [here](https://mncrowe.github.io/QGDipoles.jl/dev/Examples/)).

* Wrapper_SQG.jl: Demostrates the use of wrapper functions to easily create surface QG solutions (corresponds to Example 4 [here](https://mncrowe.github.io/QGDipoles.jl/dev/Examples/)).

### Low_Level

* 1_Layer_Modon.jl: Constructs a dipolar vortex solution in a 1-layer QG model (corresponds to Example 1 [here](https://mncrowe.github.io/QGDipoles.jl/dev/Examples/)).

* 2_Layer_Modon.jl: Constructs a dipolar vortex solution in a 2-layer QG model.

* 3_Layer_Modon.jl: Constructs a dipolar vortex solution in a 3-layer QG model (corresponds to Example 2 [here](https://mncrowe.github.io/QGDipoles.jl/dev/Examples/)).

* N_Layer_Modon.jl: Constructs a dipolar vortex solution in a general multi-layer QG model. Generalisation of 1-, 2-, and 3-layer examples.

* SQG_Modon.jl: Constructs a dipolar vortex solution in a surface QG model (corresponds to Example 3 [here](https://mncrowe.github.io/QGDipoles.jl/dev/Examples/)).

### Structures

* Vortex_Structures_LQG.jl: Constructs various 1-layer LQG vortices using vortex structures.

* Vortex_Structures_SQG.jl: Constructs various SQG vortices using vortex structures.

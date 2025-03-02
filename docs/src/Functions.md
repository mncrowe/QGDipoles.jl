# Functions

This page lists all Modules, Functions and Structures available in this package.

## QGDipoles.jl Module

```@docs
QGDipoles
```

## Vortex Structures

### LQG

```@docs
LQGParams
LQGVortex
DefLQGParams
DefLQGVortex
```

### SQG

```@docs
SQGParams
SQGVortex
DefSQGParams
DefSQGVortex
```

### Shared
```@docs
UpdateParams
UpdateVortex
```

## High-Level Functions

### LQG

```@docs
CreateModonLQG
CreateLCD
CreateLRD
```

### SQG

```@docs
CreateModonSQG
Eval_ψ_SQG
Eval_q_SQG
Eval_b_SQG
Eval_w_SQG
```

## Utilities

```@docs
CreateGrid
Calc_uv
Calc_∇
Calc_ζ
```

## Diagnostics

### LQG

```@docs
EnergyLQG
EnstrophyLQG
```

### SQG

```@docs
EnergySQG
```

## Low-Level Functions

### LQG

```@docs
BuildLinSysLQG
ApplyPassiveLayers
IncludePassiveLayers
Calc_ψq
```

### SQG

```@docs
BuildLinSysSQG
Calc_ψb
```

### Shared

```@docs
SolveInhomEVP
```

## Extras

### Monopoles

```@docs
CreateRankine
Create1LMonopole
CreateLQGMonopole
InvertVorticity1LQG
InvertVorticityLQG
```

## Internal

```@docs
QGDipoles.A_func
QGDipoles.B_func
QGDipoles.JJ_int
QGDipoles.InhomEVP_F!
QGDipoles.OrthogSpace
QGDipoles.ZernikeR
QGDipoles.GridStruct
QGDipoles.ΔNCalc
QGDipoles.CartesianGrid
QGDipoles.PolarGrid
QGDipoles.AreaInteg2
```

## Base & RecipesBase

```@docs
Base.summary
Base.show
RecipesBase.apply_recipe
```
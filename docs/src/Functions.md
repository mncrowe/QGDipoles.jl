# Functions

This page lists all Modules, Functions and Structures available in this package.

## `QGDipoles.jl`

```@docs
QGDipoles
```

## High-Level Functions

```@docs
CreateModonLQG
CreateModonSQG
CreateLCD
CreateLRD
```

## Utilities

```@docs
CreateGrid
Eval_ψ_SQG
Eval_q_SQG
Eval_b_SQG
Eval_w_SQG
Calc_∇
Calc_ζ
EnergyLQG
EnstrophyLQG
EnergySQG
```

## Vortex Types

```@docs
LQGParams
SQGParams
LQGVortex
SQGVortex
DefLQGParams
DefSQGParams
DefLQGVortex
DefSQGVortex
```

## Low-Level Functions

```@docs
BuildLinSys
ApplyPassiveLayers
IncludePassiveLayers
SolveInhomEVP
Calc_ψq
Calc_ψb
Calc_uv
```

## Monopolar Vortices (`monopoles.jl`)

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

## Base

```@docs
Base.summary
Base.show
```
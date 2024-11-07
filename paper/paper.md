---
title: 'QGDipoles.jl: A Julia package for calculating dipolar vortex solutions to the Quasi-Geostrophic equations'
tags:
  - Julia
  - rotating fluid dynamics
  - vortex dynamics
  - oceanography
  - quasi-geostrophic flows
authors:
  - name: Matthew N. Crowe
    orcid: 0000-0002-9916-2653
    equal-contrib: true
    affiliation: "1, 2"
affiliations:
 - name: School of Mathematics, Statistics and Physics, Newcastle University, Newcastle upon Tyne, NE1 7RU, UK
   index: 1
 - name: Department of Mathematics, University College London, London, WC1E 6BT, UK
   index: 2
date: 07 November 2024
bibliography: paper.bib
---

# Summary

Give summary of field, rotating flows, ocean and atmospheric modelling

# Statement of need

Use of QG flows in modelling of processes (reference), use of balanced initial conditions in full primitive equation models

# State of the field

Not much, my Matlab script, basic function from GeophysicalFlows.jl for LCD
Good recent Julia codes for QG simulations and PE simulations (GF and Oceananigans)
This package consistent with GF

# Methodology

Summarise solutions, e.g. piecewise linear relationship between (potential) vorticity and streamfunction
Give brief outline of method, emphasise scalability; i.e. coefficients calculated semi-analytically and once coefficients are found, solution can be calculated on an arbitrarily large grid (so grid independent until final stage).
Reduce differential equation to a linear algebra problem
Note current linear algebra techniques not well developed for multi-parameter problems, so root finding used

# Acknowledgements

The author would like to thank ...

# References

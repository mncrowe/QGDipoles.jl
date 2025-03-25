# Frequently Asked Questions (FAQ)

### How do I access fields inside the vortex and parameter structures?

We can access the fields of a vortex or parameter structure using a dot `.` followed by the field we want.
For example the following returns the array `b` from an `SQGVortex` structure
```julia
vortex.b
```
A full list of all fields for each structure is available on the List of Functions page.

### How do I set up solutions on a GPU?

Passing the keyword `cuda = true` to [`CreateGrid`](@ref) will create a grid on a GPU.
Any vortex solutions calculated using this grid will also be stored on the GPU.

### How do I control the accuracy of vortex solutions?

The keywords `M` and `tol` both relate to accuracy.
`M` determines the number of coefficients in the Zernike polynomial expansion used to determine the solution while `tol` determines the tolerance used when building and solving the underlying linear algebra problem.
It is recommended to set the required tolerance first, then examine the coefficients in the resulting solution.
Any coefficients less than the tolerance are set to zero so you should take the minimum `M` such that all coefficients are non-zero.
The default parameters of `tol = 1e-6` and `M = 8` (LQG) and `M = 12` (SQG) are generally sufficient.

### What grid limits can I use?

The grid limits may be set using the `Lx` and `Ly` arguments/keywords.
Setting these to scalars, e.g. `Lx = 10`, will give a grid where `x ∈ [-Lx/2, Lx/2]` and/or `y ∈ [-Ly/2, Ly/2]` and entering vector arguments, e.g. `Ly = [-4, 4]`, will give a grid where `x ∈ [Lx[1], Lx[2]]` and/or `y ∈ [Ly[1], Ly[2]]`.
We can also enter, for example, `Lx` as a scalar and `Ly` as a vector.

### How do I set the initial position and angle of the vortex?

The keywords `x₀` and `α` may be used to set the initial position and angle of the vortex.
The angle is measured anti-clockwise from the positive ``x`` axis and defined such that the vortex will travel in this direction with speed `U`.

### Can I change the method used for solving the eigenvalue problem?

Yes, you can change the method but only for single layer (including SQG) solutions.
In this case, the two available methods are
* Nonlinear root finding: root-finding is used to determine the eigenvalue and eigenvector given an initial guess.
* Generalised eigenvalue methods: the underlying problem is converted to a generalised eigenvalue problem and solved using standard methods.
The choice of method can be specified using the low-level function [`SolveInhomEVP`](@ref) with the keyword `method`.
`method = :eigensolve` corresponds to the generalised eigenvalue method while `method = :nlsolve` corresponds to nonlinear root finding.

### How do I examine the underlying linear system?

The low-level functions allow you to build the linear system explicitly.
See the [examples](@ref "Example: Low-Level Functions") here for more details.

### How do I calculate energy/enstrophy for a given vortex?

Various diagnostics can be calculated for a given vortex solution using integrals over the domain given by the `grid`.
The available diagnostics are given [here](@ref Diagnostics).
Diagnostics may be calculated directly using high-level functions or may be calculated as part of a vortex structure by setting the relevant keywords to `true` as shown [here](@ref "Example: Vortex Structures").

### How do I calculate velocity/vorticity/gradients?

The functions [`Calc_uv`](@ref), [`Calc_ζ`](@ref), and [`Calc_∇`](@ref) may be used to calculate velocities, vorticity and gradients respectively.


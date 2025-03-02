# Frequently Asked Questions (FAQ)

### How do I access fields inside the vortex and parameter structures?

We can access the fields of a vortex or parameter structure using a dot `.` followed by the field we want.
For example the following returns the array `b` from an `SQGVortex` structure
```julia
vortex.b
```
A full list of all fields for each structure is available on the List of Functions page.

### How do I control the accuracy of vortex solutions?

The keywords `M` and `tol` both relate to accuracy.
`M` determines the number of coefficients in the Zernike polynomial expansion used to determine the solution while `tol` determines the tolerance used when building and solving the underlying linear algebra problem.
It is recommended to set the required tolerance first, then examine the coefficients in the resulting solution.
Any coefficients less than the tolerance are set to zero so you should take the minimum `M` such that all coefficients are non-zero.
The default parameters of `tol = 1e-6` and `M = 8` (LQG) and `M = 12` (SQG) are generally sufficient.

### Can I change the method used for solving the eigenvalue problem?

Yes, you can change the method but only for single layer (including SQG) solutions.
In this case, the two available methods are
* Nonlinear root finding: root-finding is used to determine the eigenvalue and eigenvector given an initial guess.
* Generalised eigenvalue methods: the underlying problem is converted to a generalised eigenvalue problem and solved using standard methods.
The choice of method can be specified using the low-level function `SolveInhomEVP` with the keyword `method`.
`method = 0` corresponds to the generalised eigenvalue method while `method = 1` corresponds to nonlinear root finding.
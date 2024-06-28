# define the Zernike radial function R_n using the Jacobi polynomial

function ZernikeR(n,x)
	
	y = @. (-1)^n * jacobi(2*x^2 - 1, n, 0, 1) * (x <= 1)

	return y
end
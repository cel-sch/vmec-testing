import numpy as np
from .Bsplines import Bspline
from .library import GaussQuadrature, GaussJacobi


class Bspline_diag:
	'''
		This class creates and contains arrays with the Bsplines evaluated at the Gauss nodes for every integral. Integrals are defined in every non-zero interval.
		The dimension of the principal diagonal should be the same dimension of the matrix, which is the number of intervals plus the maximum order of the Bspline functions used:
			--> dimension = len(h) + max(nu_i,nu_j) = len(X) -1 + max(nu_i,nu_j)
			
		The class creates D matrices, with D being the total number of non-zero diagonals in the radial integration, which depends on the order of the Bspline functions used.
			--> D element [-nu_i,nu_j], D = nu_j+nu_i+1
			
		We construct 3D matrices, with dimensions (N,p,dimension-abs(diag)), with N being the number of non-zero integrals in each diagonal, p the number of Gauss nodes per integral 
		(performed on every interval) and diag the diagonal to be calculated.
			--> N = b-a, with 'a' and 'b' the integration limits on each diagonal. The interval [a,b] is of course separated in 'N' intervals on the grid.
			
		Bsplines are evaluated in a vectorised form, for an array of indexes 'alpha' or 'beta' corresponding to the Bspline index, at the corresponding 'x' values. The vectorised form
		allows for the evaluation of the full diagonal, which avoids the creation of a loop in the radial direction. We still to iterate over 3 (much smaller) loops:
			For example, assume nu_i = nu_j = 3
			Iterate over the diagonals (from -3 to 3)
				Iterate over the intervals in each diagonal (max in diag = 0: From -3 to 0)
					Iterate over the Gauss nodes (p)
					
		Note that these intervals are needed even in a non-vectorised implementation. Another possible implementation is to integrate in the full interval [a,b] with more Gauss nodes, which 
		eliminates the second loop, but might be less precise.
		
		The Bspline located in the principal diagonal appears in all the integrals. Therefore, it is evaluated in all of the intervals [-nu_i,1] and Gauss nodes. Stored in the variable 'Diag'
		
		Bsplines of Lower and Upper diagonals appear only in their respective diagonals, but spand over a few intervals. Stored in the variables Diag_L* and Diag_U*, where *=1,2,3..
		
		This loop has to be executed only once for every pair of variables nu_i,nu_j (i.e. Nvar^2 times, with Nvar the number of variables. By creating this class, we avoid the calculation of 
		these quantities Mtot^2 times, with Mtot the total number of poloidal modes.
		
		Input:
			- nu_i, nu_j : Pair of the Bspline orders
			- S          : Radial grid
			
		Attributes:
			- Diag              : Matrix containing the Bspline of the principal diagonal evaluated in the required Gauss nodes at all the intervals.
			- Diag_U*, Diag_L*  : Matrices containing the Bsplines of the Upper and Lower diagonals number (*) evaluated in the required Gauss nodes at the corresponding intervals
			
			len(Diag) = ((nu_i+1),p,dimension)
			len(Diag_U*) = len(Diag_L*) = ((b-a),p,dimension-abs(*)), where 'b' and 'a' are calculated for each diagonal (with * each diagonal number).
	
		NOTE: Some of the coefficients in the stiffness matrix have a ~1/s divergence at the axis. The Gaussian-Legendre quadrature integration rule only evaluates at the Gauss nodes, which
		do not include s=0. This means that while the numerical integration does not diverge, it is not very accurate for the very first interval where the divergence is. One could instead use the 
		Gauss-Jacobi quadrature integration rule only for the integration of the first interval, which gives a more accurate approximation of the divergent behaviour. This requires moving the nodes 
		in the first integral, which means we need to evaluate also the Bspline basis functions at the new nodes. This was done, but is only commented in the present script, the reason being that in
		several tests, a significant difference was not observed in the calculation of the eigenvalues or eigenfunctions.

	
	'''
	def __init__(self,nu_i,nu_j,der_i,der_j,grid):
		
		#Select the positiosn of the Gauss nodes.
		Gauss = GaussQuadrature(grid.NGauss)
		GaussJ = GaussJacobi(grid.NGauss)

		for diag in range(-nu_j,nu_i+1):
			
			#Bspline indexes for the rows (alpha) columns (beta) of the final matrix.
			alpha = np.arange(grid.dimension-abs(diag))-max(grid.nu)
			beta  = np.arange(grid.dimension-abs(diag))-max(grid.nu)

			a = 0
			b = min(nu_i,diag+nu_j)+1

			#Define [a,b], which are the limits of integration in each diagonal. They of course depend on the sign of diag.
			if diag < 0:
				alpha += abs(diag)
				txt = 'Diag_L'
			elif diag == 0:
				txt = 'Diag_D'
			else:
				beta += abs(diag)
				txt = 'Diag_U'
				a = diag
				

			#Create the matrices containing the diagonals
			if b>a: vars(self)[txt+str(abs(diag))] = np.zeros(((b-a),grid.NGauss,grid.dimension-abs(diag)))
		
			#Loop over the interval grid points of each diagonal (loop over the number of integrals on each diagonal)
			for l,interval in enumerate(np.arange(a,b)):
				#Calculate the value of the grid corresponding to the Bspline index (alpha) + the evaluation interval (interval).
				#For example, in an integral INT_{mu-3}^{mu-2} f(x) = 0.5*h_{mu-3}*sum_k (w_k * f(x_k))
				#with x_k = x_{mu-3}+0.5*h_{mu-3}*(1.+Pcoeff_k).
				#The variable x_ = x_{mu-3}
				#The variable h_ = h_{mu-3}
				#The variable x = x_k

				x_ = grid.S_safe[alpha+interval+grid.sk]
				h_ = grid.h_safe[alpha+interval+grid.sk]

				#Loop over the Gauss nodes
				idx0 = abs(interval+alpha[0])
				for k in range(grid.NGauss):
					x = x_ + 0.5*h_*(1.+Gauss.Pcoeff[k])
			
					# The very first interval we use Gauss-Jacobi quadrature since the integral might have a weak divergence. (uncomment if Gauss-Jacobi is desired for the first interval)
					#x[idx0] = 0.5*h_[idx0]*(1.+GaussJ.Pcoeff[k])
	
					#Calculate the product of the spline functions times the interval width. (Comment is Gauss-Jacobi is desired for the first interval)
					Result = 0.5*h_*Bspline(x,alpha,nu_i,grid,der_i)*Bspline(x,beta,nu_j,grid,der_j)*Gauss.Weight[k]

					# Uncomment if Gauss-Jacobi is desired in the first interval
					#Result = Bspline(x,alpha,nu_i,grid,der_i)*Bspline(x,beta,nu_j,grid,der_j)
					#Result[idx0] = (0.5*h_[idx0])**(1.0E-05)*Result[idx0]*GaussJ.Weight[k]
					#Result[idx0+1:] = 0.5*h_[idx0+1:]*Result[idx0+1:]*Gauss.Weight[k] 
					
					#Assing the product to the correspoinding diagonal
					vars(self)[txt+str(abs(diag))][l][k] = Result





import numpy as np
from scipy.interpolate import interp1d
import time
from scipy.fft import fftn, fftshift
from scipy.sparse import coo_matrix, diags, dok_matrix, lil_matrix
from .library import find_nearest, GaussJacobi, GaussQuadrature
from .Diagonals import Bspline_diag


def Fourier_Spline(f,hPow,grid,kind='cubic'):
	'''
		This function calculates the integral Int[ f(s,theta)*e^(i(m-m')theta) ]dtheta for defined m-m' values using a fast Fourier Transform (fft). Note that 
		one fft is performed for each radial point, so that the operation can be returned as function of 's'.
		
		Later, it performs a spline interpolation over the Gaussian nodes for the Gaussian quadrature integration in the radial direction. The results are stored
		in an array of shape = (2*Mtot-1, p, N+1), where the first dimension corresponds to the poloidal components [-(Mtot-1), (Mtot-1)].
		
		Input:
			f:		2-variable (s,theta) function to be integrated and interpolated.
			grid: 	Object containing the numerical grid in the radial and poloidal direction. It also contains several functions and arrays
					which are useful to build the matrices.
			kind:	Kind of spline interpolation (linear, quadratic, cubic, etc)
		
		Output:
			Result: Array of shape (2*Mtot-1, NGauss, N+1), containing the necessary poloidal decomposition of the function f, evaluated at the required Gauss nodes.
	'''
    
    # Fast Fourier transform on the theta component at each flux surface.
	Fourier = fftshift(fftn(f,axes=0), axes=0)/grid.Ntheta


	# Collect only the necessary Fourier modes, as we don't want to interpolate on modes that are not needed
	#-------------------------------------------------
	idx1 = find_nearest(grid.m_fft,-(grid.Mtot-1))
	idx2 = find_nearest(grid.m_fft,grid.Mtot-1)
	Fourier = Fourier[idx1:idx2+1]
	#-------------------------------------------------
	
	# Calculate the spline coefficients for all poloidal modes. Note that all the coefficients were multiplied by the function h=grid.S to eliminate the ~1/s dependence
	# in some of them. This is done because the spline interpolation does not accurately describe divergent functions. After interpolation, we divide again over s to 
	# recover the coefficients, but analytically and outside the interpolation process.
	spline = interp1d(grid.S,Fourier,kind=kind)

	# To integrate functions with weak divergences at the boundary, the Gauss-Legendre quadrature is not very precise. This is why for the very first interval one could use
	# the Gauss-Jacobi quadrature instead (commented). This doesn't change the results in any significant way (so far...) so it was decided to leave all integrations using 
	# the standard Guass-Legendre.
	GaussJ = GaussJacobi(grid.NGauss)
	Gauss = GaussQuadrature(grid.NGauss)
	
	# Interpolate poloidal modes on the required Gauss nodes for radial discretization
	# -----------------------------------------------------------------------
	Result = np.zeros(shape=(len(Fourier),grid.NGauss,grid.N+1), dtype=complex)
	for k in range(grid.NGauss):
		X_temp = grid.S[:-1] + 0.5*grid.h*(1.+Gauss.Pcoeff[k])
		
		#Comment if Gauss-Jacobi is to be used for first interval.
		Result[:,k] += spline(X_temp)*X_temp**hPow
		
		# Uncomment if Gauss-Jacobi is to be used for the first interval.
		#X_temp[0] = 0.5*grid.h[0]*(1.+GaussJ.Pcoeff[k]) 
		#Result[:,k,0] += spline(X_temp[0])
		#Result[:,k,1:] += spline(X_temp[1:])/X_temp[1:]
	# -----------------------------------------------------------------------

	return Result



def PolDisc(f,DerPol,grid,kind='cubic'):
	'''
		Let F(m,m') = Sum_{k,k'} (-i*m')^k' * (i*m)^k * Int[ f[k,k'](s,theta) * e^(i(m-m')theta) ]dtheta. This function calculates the integral of the function f[k,k'] using 
		the function 'Fourier_Spline', performs the summatory and sorts the result into a square matrix of shape=(Mtot,Mtot), where each element is multiplied by the corresponding.
		value of (-i*m')^k' * (i*m)^k.
		
		Input:
			f:			Object containing the functions f[k,k']: 2-variable (s,theta) functions to be integrated an interpolated using 'Fourier_Spline'
			DerPol:		Array containing the k,k' which are present.
			grid: 		Object containing the numerical grid in the radial and poloidal direction. It also contains several functions and arrays
						which are useful to build the matrices.
			kind:		Kind of spline interpolation (linear, quadratic, cubic, etc)
		
		Output:
			Mat: Matrix of shape=(Mtot,Mtot,p,N+1) containing all poloidal integrals of the functions f[k,k'] evaluated at the necessary Gauss nodes.
			
		NOTE: 	We concatenate zeros at the left and right of our matrix so to avoid IndexError during the radial integration process. This is harmless, as the Bsplines evaluated on
				poins outside the region [0,1] are going to be zero anyways.
	'''


	Mat_ = np.zeros(shape=(grid.Mtot,grid.Mtot,grid.NGauss,grid.N+1),dtype=complex)
	freq = np.arange(-(grid.Mtot-1),grid.Mtot)
	
	for k_ in range(len(DerPol[0])):
		# Fourier Spline integration for each f[k,k']. Note that each pair is represented by only one index 'k_'. so k = DerPol[0][k_] and k' = DerPol[1][k_]
		f_temp = Fourier_Spline(f[k_],DerPol[2][k_],grid,kind=kind)
		#f_temp = Fourier_Spline(f[k_],-1,grid,kind=kind)
		

		for i_, m_i in enumerate(grid.m):
			for j_, m_j in enumerate(grid.m):
				# Find the required poloidal mode frequency of m_i - m_j. Recall that the weight function of the fft is exp(-ikx), while our integral is exp(i(m_j-m_i)theta), and therefore the minus sign.
				idx = find_nearest(freq,m_i-m_j)
				
				if DerPol[0][k_] == 'kpar':
					'''
					In some models it's not possible to separate the terms as powers of m_i, m_j. Particularly, the parallel wave vector kpar = abs(n-m_j/q)/R.
					Here we separate the abs(ntor-m/q) term. We need to extract the following quantities from the model: 
					q    = f[k_+1][0]
					ntor = f[k_+2][0][0]
					which have been added manually to the *.py file where kpar terms are located.
					'''
					#----------------------------------------------------------
					qtemp = np.zeros(shape=(grid.NGauss,grid.N+1),dtype=complex)
					spline = interp1d(grid.S,f[k_+1][0],kind=kind)
					Gauss = GaussQuadrature(grid.NGauss)
					for l_ in range(grid.NGauss):
						X_temp = grid.S[:-1] + 0.5*grid.h*(1.+Gauss.Pcoeff[l_])
						qtemp[l_] = spline(X_temp)
					#----------------------------------------------------------
					
					Mat_[i_,j_] += (f[k_+2][0][0]+m_j/qtemp) * f_temp[idx] 
				
				else:
					Mat_[i_,j_] += (-1.0j * m_i)**DerPol[0][k_] * (1.0j * m_j)**DerPol[1][k_] * f_temp[idx]
	
	return np.concatenate([np.zeros(shape=(grid.Mtot,grid.Mtot,grid.NGauss,grid.sk)),Mat_,np.zeros(shape=(grid.Mtot,grid.Mtot,grid.NGauss,grid.sk))],axis=3)


def Discretize_Var(F_ij,grid,nu_i,nu_j,Bs_D00,Bs_D01,Bs_D10,Bs_D11,I,J,matAB,kind='cubic'):
	#Dimension of the local matrix
	localdim = grid.Mtot*grid.dimension

	#Create the full scipy matrix in coo format	
	Mat = coo_matrix((len(grid.nu)*localdim,len(grid.nu)*localdim))
	
	if F_ij.DerBs[0]: F_D00 = PolDisc(F_ij.var_D00, F_ij.DerPol_DerBs00, grid, kind=kind)
	if F_ij.DerBs[1]: F_D01 = PolDisc(F_ij.var_D01, F_ij.DerPol_DerBs01, grid, kind=kind)
	if F_ij.DerBs[2]: F_D10 = PolDisc(F_ij.var_D10, F_ij.DerPol_DerBs10, grid, kind=kind)
	if F_ij.DerBs[3]: F_D11 = PolDisc(F_ij.var_D11, F_ij.DerPol_DerBs11, grid, kind=kind)
	
	#Loop over the diagonal elements
	for diag in range(-nu_j,nu_i+1):

		#Indexes where to evaluate the functions
		alpha = np.arange(grid.dimension-abs(diag),dtype=np.int32)-max(grid.nu)
		a = 0
		b = min(nu_i,diag+nu_j)+1
		
		if diag < 0:
			alpha += abs(diag)
			txtD = 'Diag_L'
			txti = 'I_L'
			txtj = 'J_L'
			
		elif diag == 0:
			txtD = 'Diag_D'
			txti = 'I_D'
			txtj = 'J_D'
		else:
			a = diag
			txtD = 'Diag_U'
			txti = 'I_U'
			txtj = 'J_U'
			
		#Each integral is defined in an Interval [a,b]. Each interval can be composed of many sub-intervals, corresponding to the grid intervals.
		#Example: [a,b] = [-1,3] = [-1,0] + [0,1] + [1,2] + [2,3]
		#The following loop is over the sub-intervals.
		DiagsAll = np.zeros(shape=(grid.Mtot,grid.Mtot,grid.dimension-abs(diag)),dtype=complex)
		
		for l,interval in enumerate(np.arange(a,b)):
			#Evaluate the functions in the required interval
			if F_ij.DerBs[0]: F_D00_ = F_D00[:,:,:,alpha+interval+grid.sk]
			if F_ij.DerBs[1]: F_D01_ = F_D01[:,:,:,alpha+interval+grid.sk]
			if F_ij.DerBs[2]: F_D10_ = F_D10[:,:,:,alpha+interval+grid.sk]
			if F_ij.DerBs[3]: F_D11_ = F_D11[:,:,:,alpha+interval+grid.sk]

			#Loop over the Gaussian nodes at each sub-interval
			for k in range(grid.NGauss):

				if F_ij.DerBs[0]: DiagsAll += vars(Bs_D00)[txtD+str(abs(diag))][l][k]*F_D00_[:,:,k]
				if F_ij.DerBs[1]: DiagsAll += vars(Bs_D01)[txtD+str(abs(diag))][l][k]*F_D01_[:,:,k]
				if F_ij.DerBs[2]: DiagsAll += vars(Bs_D10)[txtD+str(abs(diag))][l][k]*F_D10_[:,:,k]
				if F_ij.DerBs[3]: DiagsAll += vars(Bs_D11)[txtD+str(abs(diag))][l][k]*F_D11_[:,:,k]
					

		#Account for mixed finite element basis
		#-------------------------------------------------------------------------------
		#if matAB == 'B':
		if I == J and diag == 0 and (nu_i < max(grid.nu) or nu_j < max(grid.nu)):
			for mix in range(max( abs(nu_i-max(grid.nu)), abs(nu_j-max(grid.nu)))):
				for mmix in range(grid.Mtot):
					DiagsAll[mmix,mmix,mix] = 1.
		#-------------------------------------------------------------------------------
		
		#Set the values directly on the PETSc matrix using the RowColmnValue format
		#--------------------------------------------------------------------------
		Mat += coo_matrix((np.reshape(DiagsAll,grid.Mtot**2*(grid.dimension-abs(diag))),(localdim*I+vars(grid)[txti+str(abs(diag))],localdim*J+vars(grid)[txtj+str(abs(diag))])),shape=(len(grid.nu)*localdim,len(grid.nu)*localdim))
		#--------------------------------------------------------------------------
		
	return Mat





#=========================================================================================================================================
def Discretize_Operators_scipy(grid,EQUIL,ntor,Coeff_A,Coeff_B,Write=False):
	'''
		Global function to discretize the operators A and B, which constitute the general eigenvalue problem A*v = w B*v. 
		
		Input:
			grid: 		Object containing the numerical grid in the radial and poloidal direction. It also contains several functions and arrays
						which are useful to build the matrices.
			EQUIL:		Structure containing all the necessary equilibrium quantities in SFL coordinates.
			Coeff_A,B:	Objects (stored in a separate file usually) containing the operators to be applied to the vector v. This defines the model
	'''

	#Dimension of the matrix of each variable.
	localdim = grid.Mtot*grid.dimension

	# Create and allocate matrices A and B
	t0 = time.time()
	#-----------------------------------------------------------------------------------------------
	A = coo_matrix((len(grid.nu)*localdim,len(grid.nu)*localdim))
	B = coo_matrix((len(grid.nu)*localdim,len(grid.nu)*localdim))
	#-----------------------------------------------------------------------------------------------
	
	# We calculate the Bspline_diag only for different nu_i,nu_j pairs, then we discretise all the variables with equal nu_i,nu_j pair.
	# This avoids repeating the calculation of Bspline_diag for equal nu_i,nu_j pairs.
	#----------------------------------------------------------------------------------------------------------------------------------
	# Loop over (nu_i,nu_j) pairs
	for i_pair,pair in enumerate(grid.pair):
		
		nu_i = pair[0]
		nu_j = pair[1]
		
		# Calculate the Bspline elements for each combination of derivatives (considering 1 to be the highes order derivative in the weak form)
		# Note that not all elements are needed in a given weak formulation. I'll leave them here for generality, but one can comment that lines 
		# that are not needed to make the code a little bit faster.
		Bs_D00 = Bspline_diag(nu_i,nu_j,0,0,grid)
		Bs_D01 = Bspline_diag(nu_i,nu_j,0,1,grid)
		Bs_D10 = Bspline_diag(nu_i,nu_j,1,0,grid)
		Bs_D11 = Bspline_diag(nu_i,nu_j,1,1,grid)

		
		# Loop over the variables using the same (nu_i, nu_j) pairs. 
		for loc in grid.indices[i_pair]:
			
			#Indices connecting the submatrices with the global matrix.
			i = loc[0]
			j = loc[1]
			
			#Calculate the coefficients related to this variable.
			A_ij  = vars(Coeff_A)['M_'+str(i)+str(j)](EQUIL,ntor)
			B_ij  = vars(Coeff_B)['M_'+str(i)+str(j)](EQUIL,ntor)
			
			# Discretize operators A and B
			#--------------------------------------------------------------------------------------------------------------
			A += Discretize_Var(A_ij,grid,nu_i,nu_j,Bs_D00,Bs_D01,Bs_D10,Bs_D11,i,j,'A')
			B += Discretize_Var(B_ij,grid,nu_i,nu_j,Bs_D00,Bs_D01,Bs_D10,Bs_D11,i,j,'B')
			#--------------------------------------------------------------------------------------------------------------
			

	print ('Discretisation finished in %.5E s'%(time.time()-t0))

	#Apply boundary conditions
	#3 Types of BC can be applied: Dirichlet, Neumann and Homogeneous Robin (mixed). Homogeneous Robin BC are of the form alpha*dX/ds+beta*X = 0,
	#which one can write as a non-homogeneous Neumann BC dX/ds = -beta/alpha * X, and apply it as a Natural BC using the boundary terms obtained
	#by integrations by parts in the proces of obtaining the weak form. Homogeneous Robin BC can be of the form alpha*dX/ds+beta*Y=0, where X and 
	#Y are different variables. In that case, one can substitute dX/ds=-beta/alpha * Y but now the terms are crossed: BC needs to be applied to 
	#variable Y, but in the equation of variable X. Therefore, for Robin BC one must specify with which variables the BC are mixed (can be more 
	#than one!, so careful). BC are specified in the file with the coefficients of the matrix A.
	#----------------------------------------------------------------------------------------------------------------------------------
	t0 = time.time()
	
	#Transform to lil format to easy access to the elements
	A = A.tolil()
	B = B.tolil()
	
	Dirichlet_idx = []
	Neumann_idx = []

	for i_bc, nu in enumerate(grid.nu):
		# Location of the boundary terms
		x0 = np.arange(i_bc*localdim,(i_bc+1)*localdim, grid.dimension, dtype=int)+abs(max(grid.nu)-nu)
		x1 = x0+grid.dimension-1-abs(max(grid.nu)-nu)

		BC = vars(Coeff_A)['BC_'+str(i_bc)](grid,EQUIL,ntor)

		for j_bc in range(grid.Mtot):
			#BC at the magnetic axis
			if BC.Axis[j_bc] == 'Dirichlet':
				Dirichlet_idx.append(x0[j_bc])

			elif BC.Axis[j_bc] == 'Neumann':
				Neumann_idx.append(x0[j_bc])

			elif BC.Axis[j_bc] == 'Robin':
				k_bc = BC.RobinAxisMixed
				x0_ = np.arange(k_bc*localdim,(k_bc+1)*localdim, grid.dimension, dtype=int)+abs(max(grid.nu)-nu)
				if i_bc <= k_bc: A[x0[j_bc],x0_[j_bc]] += BC.RobinAxisValue[j_bc]
				if i_bc >  k_bc: A[x0_[j_bc],x0[j_bc]] += BC.RobinAxisValue[j_bc]
			else:
				None
			
			#BC at the plasma edge
			if BC.Edge[j_bc] == 'Dirichlet':
				Dirichlet_idx.append(x1[j_bc])
				
			elif BC.Edge[j_bc] == 'Neumann':
				Neumann_idx.append(x1[j_bc])
				
			elif BC.Edge[j_bc] == 'Robin':
				k_bc = BC.RobinEdgeMixed
				x1_ = np.arange(k_bc*localdim,(k_bc+1)*localdim, grid.dimension, dtype=int)+grid.dimension-1
				if i_bc <= k_bc: A[x1[j_bc],x1_[j_bc]] += BC.RobinEdgeValue[j_bc]
				if i_bc >  k_bc: A[x1_[j_bc],x1[j_bc]] += BC.RobinEdgeValue[j_bc]
			else:
				None
	
	#Apply the boundary conditions
	#--------------------------------------------------
	Dirichlet_idx = np.asarray(Dirichlet_idx,dtype=int)
	Neumann_idx = np.asarray(Neumann_idx,dtype=int)
	A = Dirichlet(A,Dirichlet_idx,1.)
	B = Dirichlet(B,Dirichlet_idx,0.)

	A = Neumann(A,Neumann_idx,1.)
	B = Neumann(B,Neumann_idx,0.)
	#--------------------------------------------------

	print ('Time to apply boundary conditions: %.5E s'%(time.time()-t0))
	
	#----------------------------------------------------------------------------------------------------------------------------------

	return A.tocsc(),B.tocsc()
#=========================================================================================================================================



def Dirichlet(mat,idx,value_diag=1.):

	mat[idx,:] = 0.
	mat[:,idx] = 0.
	mat[idx,idx] = value_diag

	return mat

def Neumann(mat,idx,value_diag=1.):
    
	mat[idx,:] = 0.
	mat[idx,idx] = value_diag
	mat[idx,idx+1] = -value_diag
		
	return mat




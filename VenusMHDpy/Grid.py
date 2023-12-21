import numpy as np
from .library import find_nearest

class GRID:
	'''
		This class should contain all the quantities related to the numerical grid used in the discretization process. The parameters needed to generate the full grid are:
			
			N          : The number of interior grid points. Total number of points should be N+2.
			nu         : An array of the Bspline orders of the different variables. This does not modify the grid, but defines the size of the matrices, how sparse they are and the structure of
						 the global matrix.
			Mmin, Mmax : Minimum and Maximum poloidal modes used in the Fourier discretization in the poloidal direction.
			NGauss     : Number of Gauss nodes to be used in the radial integration process.
			bunching   : Whether increase the density of radial points around certain plasma radii (at rational surfaces for example)
	'''
	def __init__(self):
		#Default values
		
		#Radial grid
        #---------------------
		self.N = 500
		self.nu = [3,2,2]
		self.S = np.linspace(0.,1.,self.N+2)
		
		#Poloidal grid
        #---------------------
		self.Mmin = -3
		self.Mmax = 5
		self.Ntheta = 50

		#Number of Gauss Nodes
		#---------------------
		self.NGauss = 5
		
		#Bunching variables
		#---------------------
		self.bunching = False
		self.bunchingValues = [0.5]
		self.bunchingQValues =[1.]
		self.bunchingAmplitudes = [2.0]
		self.bunchingSigma = [0.01]
		self.qtol = 1.0E-02

	def BuildGrid(self, s_qfactor=None, qfactor=[]):
		'''
			Builds the grid and the global matrix  structure.
			If bunching = True, then we must specify either:
				-The q factor, the associated radial coordinate and the rational surfaces where bunching will occur.
				-The radial points where bunching will occur.
					-Additionally, one must specify the amplitudes and sigma values for each bunching point.
		'''
		#Radial grid
		#---------------------------------------------------
		if self.bunching:
			if len(qfactor) == 0:
				
				self.S = self.CreateBunching(self.N,self.bunchingAmplitudes,self.bunchingValues,self.bunchingSigma)
				self.N = len(self.S)-2
			
			else:
				bunchingValues = []
				bunchingAmplitudes = []
				bunchingSigma = []
	
				for l,q_ in enumerate(self.bunchingQValues):
					idx_ = find_nearest(qfactor,q_)
					if abs(qfactor[idx_]-q_) < self.qtol:
						bunchingValues.append(s_qfactor[idx_])
						bunchingAmplitudes.append(self.bunchingAmplitudes[l])
						bunchingSigma.append(self.bunchingSigma[l])
					else:
						print ('Rational surface %.4f not found in safety factor for bunching'%(self.bunchingQValues[l]))

				if len(bunchingValues) > 0:
					self.S = self.CreateBunching(self.N,bunchingAmplitudes,bunchingValues,bunchingSigma)
					self.N = len(self.S)-2
				else:
					self.S = np.linspace(0.,1.,self.N+2)

		else:
			#delta = 1.0E+06
			self.S = np.linspace(0.,1.,self.N+2)
		
		self.dimension = self.N+1+max(self.nu)
		self.h = np.diff(self.S)
		
		# To avoid calculating off-grid operations separately from the bulk, it is useful to add dummy values at the beggining and at the end of the arrays. This way, one can perform all the multiplications simultaneusly
		# When indexing, the indices will be shifted by the number of dummy values added (sk). A safe number of sk is max(nu) + 1.
		self.sk = max(self.nu)+1
		self.knots  = np.concatenate([np.ones(self.sk,dtype=np.int64)*self.S[0],self.S,np.ones(self.sk,dtype=np.int64)*self.S[-1]])
		self.S_safe = np.concatenate([np.zeros(self.sk,dtype=np.int64),self.S,np.zeros(self.sk,dtype=np.int64)])
		self.h_safe = np.concatenate([np.zeros(self.sk,dtype=np.int64),self.h,np.zeros(self.sk,dtype=np.int64)])
		
		#Poloidal grid
		#---------------------------------------------------
		#One should be capable to resolve freq.uencies from -(Mtot-1) to +(Mtot-1). So in theory one should have Ntheta = 2*(Mtot-1)+1. Just in case, we set Ntheta = Mtot*4+1.
		self.Mtot = self.Mmax-self.Mmin+1
		# self.Ntheta = self.Mtot*4+1
		self.theta = np.linspace(0.,2.*np.pi,self.Ntheta,endpoint=False)
		self.m = np.arange(self.Mmin,self.Mmax+1)
		self.m_fft = np.arange( -int(np.fix(self.Ntheta/2)),int(np.ceil(self.Ntheta/2)))
		
		
		#---------------------------------------------------
		self.NNZ = 0
		for diag_ in range(-max(self.nu),max(self.nu)+1):
			'''
				Set arrays with coo indices with diagonal elements. This way one can set each diagonal for all radial points in all poloidal modes for each variable using only one assigment command.
				CSR might be faster if we use local matrices with PETSc BAIJ matrix. BUT the block matrices are apparently dense by default, and I haven't been able to modify this.
			'''
			I,J = self.coo_idx(diag_)
			
			if diag_ < 0: 
				vars(self)['I_L'+str(abs(diag_))] = I
				vars(self)['J_L'+str(abs(diag_))] = J
			if diag_ > 0: 
				vars(self)['I_U'+str(abs(diag_))] = I
				vars(self)['J_U'+str(abs(diag_))] = J
			if diag_ == 0: 
				vars(self)['I_D'+str(abs(diag_))] = I
				vars(self)['J_D'+str(abs(diag_))] = J
				
		self.MatStruct_nu()
				
    
	def MatStruct_nu(self):
		'''
			Map of the global matrix structure, indicating which submatrices use the same Bspline combination. 
			Works as expected but not very elegant. Anyways, this has to be run only once in the whole code.
		'''
		self.pair = []
		self.indices = []
		for i,nu_i in enumerate(self.nu):
			for j,nu_j in enumerate(self.nu):
				if [nu_i,nu_j] in self.pair:
					for k,v in enumerate(self.pair):
						if [nu_i,nu_j] == v:
							self.indices[k].append([i,j])
							
				else:
					self.pair.append([nu_i,nu_j])
					self.indices.append([])
					self.indices[-1].append([i,j])


	def coo_idx(self,diag):
		'''
			Calcuates the local coo indexes (I,J) for a given diagonal (diag)
		
		'''

		#Indexes for the principal diagonal (alpha) and upper and lower diagonals (beta).
		alpha = np.arange(self.dimension-abs(diag), dtype=np.int64)
		beta  = np.arange(self.dimension-abs(diag), dtype=np.int64)

		if diag < 0:
			alpha += abs(diag)
		else:
			beta += abs(diag)
			
		I = []
		J = []
		for m_ in range(self.Mtot):
			I = np.append(I,np.tile(alpha+m_*self.dimension,self.Mtot))
			J = np.append(J,beta+m_*self.dimension)
			
		J = np.tile(J,self.Mtot)

		# return np.reshape(I.astype(np.int64),(len(I),1)),np.reshape(J.astype(np.int64),(len(J),1))
		return I.astype(np.int64),J.astype(np.int64)
	
	def getNNZ_local(self,nu_i,nu_j):
		'''
			Calculates the number of elements in a local matrix whose variables are expressed in Bsplines of orders (nu_i,nu_j).
			This is to reserve the amount of memory needed to build a local matrix, which is useful if PETSc is used.

		'''
		return int((self.dimension*(nu_i+nu_j+1)-nu_i*(nu_i+1)/2-nu_j*(nu_j+1)/2)*self.Mtot**2)

	def getNNZ_global(self):
		'''
			Calculates the number of elements in the full matrix. Useful if PETSc is used.
		'''
		nnz = 0
		for i in self.nu:
			for j in self.nu:
				nnz += self.getNNZ_local(i,j)
		return nnz
    
	def BunchingDensity(self,x,A,r0,sigma):
		"""
			Function defining the density of points. The density function is a Gaussian where:
				A    : Amplitude of the Gaussian
				r0   : Location of the Gaussian
				sigma: Standard deviation of the Gaussian
			
			Note that A,r0 and sigma must be arrays with the same number of elements.
		"""
		
		fdens = 1.
		for i in range(len(r0)):
			
			fdens += A[i]*np.exp(-0.5*((x-r0[i])/sigma[i])**2.)
		
		return fdens
	
	def CreateBunching(self,N,A,r0,sigma, ds_init=0.01):
		"""
			Creates the radial grid with bunching around certain points.
		"""
		
		x_ = np.linspace(0.,1.,10000)
		fdens = self.BunchingDensity(x_,A,r0,sigma)
		
		fdens0 = N/np.trapz(fdens,x_)
		
		f = [0.]
		ds = ds_init
		while f[-1] < 1.:
			k = 0.5*fdens0*(self.BunchingDensity(f[-1],A,r0,sigma)+self.BunchingDensity(f[-1]+ds,A,r0,sigma))*ds
			ds = ds/k
			
			if f[-1]+ds < 1.0: 
				f.append(f[-1]+ds)
			else:
				if 1.-f[-1] < f[-1]+ds-1.:
					f[-1] = 1.
				else:
					f.append(1.0)
	
		return np.asarray(f)






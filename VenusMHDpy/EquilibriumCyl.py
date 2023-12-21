import numpy as np 
from .library import GaussQuadrature

'''
	Cylindrical Equilibria

	Equilibrium quantities for plasmas in cylindrical geometry. Used to benchmark the code. The examples listed below are analytical equilibria obtained from the book
	"Finite Element Methods in Linear Ideal Magnetohydrodynamics" by Ralf Gruber and Jacques Rappaz.

	Last modified: 08/April/2022 
'''


class Profiles:
	def __init__(self,grid,jz=0.):		
		
		self.R = np.ones(shape=(grid.Ntheta,len(grid.S)))
		self.gamma = 5./3.

		Gauss = GaussQuadrature(grid.NGauss)
		self.S = grid.S.copy()
		
		
		# Test case F
		#-------------------------------------------------------------------------------
		c1 = 2./7.
		c2 = 10./7.
		self.Bt = c1*self.S/(1.+c2**2.*self.S**2.)
		self.dBtds = c1*(1.-c2**2.*self.S**2.)/(1.+c2**2.*self.S**2.)**2.
		self.Bz = np.ones_like(self.S)
		self.dBzds = np.zeros_like(self.S)
		
		self.rho0 = np.ones_like(self.S)
		self.P = 0.5*(c1/c2)**2.*(1./(1.+c2**2.*self.S**2.)**2. - 1./(1.+c2**2.)**2.)
		
		self.q = (1.+c2**2.*self.S**2.)/c1
		self.dqds = 2.*self.S*c2**2./c1
		#-------------------------------------------------------------------------------
				
		# Test case D
		#-----------------------------------------
		#self.jz = jz
		#self.Bt = 0.5*grid.S*self.jz
		#self.dBtds = 0.5*self.jz*np.ones_like(grid.S)
		#self.Bz = np.ones_like(grid.S)
		#self.dBzds = np.zeros_like(grid.S)
		
		#self.rho0 = np.ones_like(grid.S)
		#self.P = 0.25*self.jz**2.*(1.-grid.S**2.)
		
		#self.q = 2.*np.ones_like(grid.S)/self.jz
		#self.dqds = 0.
		#-----------------------------------------

		
				
		

		
		
		
		

		


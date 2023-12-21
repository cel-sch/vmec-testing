import numpy as np
import sys
sys.path.insert(0, '../bin/')
import library as lib



class M_00:
	def __init__(self,grid,ntor):
		
		#Convention is: 00,01,10,11.
		self.DerBs = [True,False,False,False]
		
		self.DerPol_DerBs00 = [[0],[0]]
		
		self.var_D00 = np.zeros(shape=(len(self.DerPol_DerBs00[0]),grid.Ntheta, grid.N+2))
		
		#First set of coefficients
		#self.var_D00[0] = grid.S**3.

		#Second set of coefficients
		self.var_D00[0] = grid.S**2.

import numpy as np
import sys
sys.path.insert(0, '../bin/')
import library as lib



class M_00:
	def __init__(self,grid,ntor):
		
		#Convention is: 00,01,10,11.
		self.DerBs = [True,True,False,True]
		#self.DerBs = [True,False,False,True]
		
		# self.DerPol_DerBs00 = [[0],[0]]
		self.DerPol_DerBs00 = [[0],[2]]
		self.DerPol_DerBs01 = [[0],[0]]
		self.DerPol_DerBs11 = [[0],[0]]
		
		self.var_D00 = np.zeros(shape=(len(self.DerPol_DerBs00[0]),grid.Ntheta, grid.N+2))
		self.var_D01 = np.zeros(shape=(len(self.DerPol_DerBs01[0]),grid.Ntheta, grid.N+2))
		self.var_D11 = np.zeros(shape=(len(self.DerPol_DerBs11[0]),grid.Ntheta, grid.N+2))

	
		#First set of coefficients
		#self.var_D00[0] = -grid.S
		#self.var_D01[0] = grid.S**2.
		#self.var_D11[0] = grid.S**3.
		
		#Second set of coefficients
		self.var_D00[0] = -np.ones_like(grid.S)
		self.var_D11[0] = grid.S**2.
		
		#NOTE 
		#The first set of coefficients correspond to solving Int( z* (s d2wds + dwds + dwdu/s = -lambda s w) ) s dsdu. This representation doesn't give 1/s coefficients.
		#The second set of coefficients correspond to solving Int( z* (d2wds + dwds/s + dwdu/s**2 = -lambda w) ) s dsdu. This representation gives 1/s coefficients.
		#Convergence order for both representations is 2*nu.
		
class BC_0:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Dirichlet']*grid.Mtot)
       
		#Normally the m=0 mode should have natural BC 
        idx = np.where(abs(grid.m) == 0)[0]
        self.Axis[idx] = ['Natural']

import numpy as np
import sys
sys.path.insert(0, '../bin/')




class M_00:
	def __init__(self,grid,ntor):
		
		#Convention is: 00,01,10,11.
		self.DerBs = [True,False,False,True]
		
		self.DerPol_DerBs00 = [[0],[0]]
		self.DerPol_DerBs11 = [[0],[0]]
		
		self.var_D00 = np.zeros(shape=(len(self.DerPol_DerBs00[0]),grid.Ntheta, grid.N+2))
		self.var_D11 = np.zeros(shape=(len(self.DerPol_DerBs11[0]),grid.Ntheta, grid.N+2))
		
		# alpha = 1.+grid.S
		alpha = 1.
		gamma = np.zeros_like(grid.S)
		
		self.var_D00[0] = gamma*grid.S
		self.var_D11[0] = alpha*grid.S
		
			
class M_01:
	def __init__(self,grid,ntor):
		
		#Convention is: 00,01,10,11.
		self.DerBs = [True,True,False,False]
		
		self.DerPol_DerBs00 = [[0],[0]]
		self.DerPol_DerBs01 = [[0],[0]]
		
		self.var_D00 = np.zeros(shape=(len(self.DerPol_DerBs00[0]),grid.Ntheta, grid.N+2))
		self.var_D01 = np.zeros(shape=(len(self.DerPol_DerBs01[0]),grid.Ntheta, grid.N+2))
		
		beta = np.ones_like(grid.S)
		betap = np.zeros_like(grid.S)
		
		self.var_D00[0] = -betap*grid.S
		self.var_D01[0] = -beta*grid.S

		
class M_10:
	def __init__(self,grid,ntor):
		
		#Convention is: 00,01,10,11.
		self.DerBs = [False,True,False,False]
		
		self.DerPol_DerBs01 = [[0],[0]]
		
		self.var_D01 = np.zeros(shape=(len(self.DerPol_DerBs01[0]),grid.Ntheta, grid.N+2))
		
		beta = np.ones_like(grid.S)

		self.var_D01[0] = beta*grid.S
			
class M_11:
	def __init__(self,grid,ntor):
		
		#Convention is: 00,01,10,11.
		self.DerBs = [True,False,False,False]
		
		self.DerPol_DerBs00 = [[0],[0]]
		
		self.var_D00 = np.zeros(shape=(len(self.DerPol_DerBs00[0]),grid.Ntheta, grid.N+2))
		
		delta = np.ones_like(grid.S)
		
		self.var_D00[0] = delta*grid.S


# Boundary conditions
#==============================================================#


class BC_0:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Robin']*grid.Mtot)

        self.RobinEdgeMixed = 1
        self.RobinEdgeValue = np.ones(grid.Mtot)
        
class BC_1:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)

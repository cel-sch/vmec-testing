import numpy as np
import sys
sys.path.insert(0, '../bin/')


class M_00:
    def __init__(self,grid,ntor):
        
        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]
        
        self.DerPol_DerBs00 = [[0],[0]]
        
        self.var_D00 = np.zeros(shape=(len(self.DerPol_DerBs00[0]),grid.Ntheta, grid.N+2))
        
        rho = np.ones_like(grid.S)
        
        self.var_D00[0] = rho*grid.S
        
class M_01:
    def __init__(self,grid,ntor):
        
        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]
	
            
class M_10:
    def __init__(self,grid,ntor):
        
        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]


class M_11:
    def __init__(self,grid,ntor):
        
        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]
        
        self.DerPol_DerBs00 = [[0],[0]]
        
        self.var_D00 = np.zeros(shape=(len(self.DerPol_DerBs00[0]),grid.Ntheta, grid.N+2))
        
        rho = np.ones_like(grid.S)
        
        self.var_D00[0] = rho*grid.S

import numpy as np 
from .library import GaussQuadrature

'''
	Parameters for the infinite-n ballooning equation.
'''


class Profiles:
	def __init__(self,grid):		
		
		self.R = np.ones(shape=(grid.Ntheta,len(grid.S)))

		self.r = grid.S.copy()
		self.alpha = 1.
		self.shear = 1.
		self.dm = -0.2
		self.k = 0.
		self.delta = 1.0E-09

import numpy as np 
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev, interp1d
from .library import SplineDer1, SplineDer2

class Profiles:
	def __init__(self,SFLFile):
		File = h5py.File(SFLFile, 'r')
		#Recall f=g*h, where h=grid.S

		self.s = np.asarray(list(File['profiles']['s']))
		self.q = np.asarray(list(File['profiles']['q']))
		self.T = np.asarray(list(File['profiles']['T']))
		self.F = np.asarray(list(File['profiles']['F']))
		self.g = np.asarray(list(File['profiles']['g']))
		self.h = np.asarray(list(File['profiles']['h']))
		self.P = np.asarray(list(File['profiles']['P']))
		self.R = np.asarray(list(File['geometry']['R']))
		self.Z = np.asarray(list(File['geometry']['Z']))
		self.R0 = File['normalisation']['R0'][()]
		self.B0 = File['normalisation']['B0'][()]
		self.rho0 = np.ones_like(self.h)				#Constant density (Normalised)
		self.rho = np.ones_like(self.h)				#Constant density (Normalised)
		#self.rho0 = self.P/self.P[0]					#Constant temperature (Normalised)

		self.gamma = 5./3. 
		self.mu0 = 4.*np.pi*1.0E-07

	def Normalise(self):
		
		self.F = self.F/(self.R0*self.B0)
		self.g = self.g/(self.R0**2.*self.B0)
		self.P = self.mu0*self.P/self.B0**2.
		self.T = self.mu0*self.T/self.B0**2.
		self.R = self.R/self.R0
		self.Z = self.Z/self.R0
		

	def ChangeGrid(self,S_new):
	
		dummy  = splrep(self.s,self.q)
		self.q = splev(S_new,dummy)

		dummy  = splrep(self.s,self.F)
		self.F = splev(S_new,dummy)
		
		dummy  = splrep(self.s,self.g)
		self.g = splev(S_new,dummy)
		
		dummy  = splrep(self.s,self.h)
		self.h = splev(S_new,dummy)

		dummy  = splrep(self.s,self.P)
		self.P = splev(S_new,dummy)

		dummy  = splrep(self.s,self.rho0)
		self.rho0 = splev(S_new,dummy)

		dummy  = interp1d(self.s,self.R,kind='cubic',axis=1)
		self.R = dummy(S_new)   

		dummy  = interp1d(self.s,self.Z,kind='cubic',axis=1)
		self.Z = dummy(S_new)
		
		self.s = S_new

	def BuildInGrid(self,grid):
	
		#Variables
		#===================================
		self.dgds   = SplineDer1(self.g,grid,der_s=1)
		self.dhds   = SplineDer1(self.h,grid,der_s=1)
		self.dqds   = SplineDer1(self.q,grid,der_s=1)
		self.d2qds  = SplineDer1(self.q,grid,der_s=2)
		self.dPds   = SplineDer1(self.P,grid,der_s=1)
		self.d2Pds  = SplineDer1(self.P,grid,der_s=2) 
		self.dFds   = SplineDer1(self.F,grid,der_s=1)
		self.dRds   = SplineDer2(self.R,grid,der_s=1,der_theta=0)
		self.dRdu   = SplineDer2(self.R,grid,der_s=0,der_theta=1)
		self.dZds   = SplineDer2(self.Z,grid,der_s=1,der_theta=0)
		self.dZdu   = SplineDer2(self.Z,grid,der_s=0,der_theta=1)
		self.R2     = self.R*self.R
		self.dR2ds  = SplineDer2(self.R2,grid,der_s=1,der_theta=0)
		self.dR2du  = SplineDer2(self.R2,grid,der_s=0,der_theta=1)
		self.g11    = self.dRds**2 + self.dZds**2
		self.dg11ds = SplineDer2(self.g11,grid,der_s=1,der_theta=0)
		self.dg11du = SplineDer2(self.g11,grid,der_s=0,der_theta=1)
		self.g12    = self.dRds*self.dRdu + self.dZds*self.dZdu
		self.dg12ds = SplineDer2(self.g12,grid,der_s=1,der_theta=0)
		self.dg12du = SplineDer2(self.g12,grid,der_s=0,der_theta=1)
		self.g22    = self.dRdu**2 + self.dZdu**2
		self.dg22ds = SplineDer2(self.g22,grid,der_s=1,der_theta=0)
		self.dg22du = SplineDer2(self.g22,grid,der_s=0,der_theta=1)
		#self.J2 = -((self.dFds*self.F)/(self.g*self.h*self.q*self.R2))
		#self.J3 = -1./(self.g*self.h)*(self.F*self.dFds/self.R2+self.dPds)
		
		self.B2     = self.F**2.*(self.q**2.*self.R2+self.g22)/(self.q**2.*self.R2**2)
		self.dB2ds  = SplineDer2(self.B2,grid,der_s=1,der_theta=0)
		self.dB2du  = SplineDer2(self.B2,grid,der_s=0,der_theta=1)
		
		imag = 1.0j
		#===================================


	def plot(self,grid,show=True):

		ax1 = plt.subplot(231)
		plt.plot(grid.S,self.F*self.R0*self.B0)
		plt.tick_params('x', labelbottom=False)
		plt.ylabel(r'F $[Wb/m]$')
		plt.xlabel('s')
		plt.grid()

		ax2 = plt.subplot(232)
		plt.plot(grid.S,self.g*self.h*self.R0**2.*self.B0)
		plt.tick_params('x', labelbottom=False)
		plt.xlabel('s')
		plt.ylabel('f [Wb]')
		plt.grid()

		ax3 = plt.subplot(234)
		plt.plot(grid.S,self.B0**2.*self.P/self.mu0)
		plt.tick_params('x', labelbottom=False)
		plt.xlabel('s')
		plt.ylabel('P [Pa]')
		plt.grid()
	
		ax5 = plt.subplot(235)
		plt.plot(grid.S,self.q)
		plt.ylabel(r'q')
		plt.xlabel('s')
		plt.grid()

		ax6 = plt.subplot(133)
		plt.plot(self.R[:,-1]*self.R0,self.Z[:,-1]*self.R0)
		plt.plot(self.R[0,0]*self.R0,self.Z[0,0]*self.R0,'+')
		plt.xlabel(r'R [m]')
		plt.ylabel(r'Z [m]')
		plt.gca().set_aspect('equal')
		plt.grid()
		
		
		#Check Grad Shafranov
		#LHS = (SplineDer2(self.F*self.g22/(self.q*self.R2),grid,der_s=1,der_theta=0)-SplineDer2(self.F*self.g12/(self.q*self.R2),grid,der_s=0,der_theta=1))
		LHS = self.dFds*self.g22/(self.q*self.R2)+self.F*self.dg22ds/(self.q*self.R2)-self.F*self.g22*self.dqds/(self.q**2*self.R2)-self.F*self.g22*self.dR2ds/(self.q*self.R2**2.)-self.F/self.q*(self.dg12du/self.R2-self.g12*self.dR2du/self.R2**2.)
		RHS = -self.q*self.R2/self.F*(self.F*self.dFds/self.R2+self.dPds)
		
		x = 0
		plt.figure(2)
		plt.plot(grid.S,LHS[x], 'x-', label=r'$\Delta^*\psi_p/R^2$')
		plt.plot(grid.S,RHS[x], 'o-', label=r'$J^3$')
		plt.plot(grid.S,(RHS[x]-LHS[x])/max(abs(LHS[x])), 'o-', label=r'$\frac{|J^3-\Delta^*\psi_p/R^2|}{max(\Delta^*\psi_p/R^2)}$')
		plt.grid()

		plt.legend()
		if show: plt.show()



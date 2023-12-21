import numpy as np 
import matplotlib.pyplot as plt

class Profiles:
	def __init__(self,R0,B0,q0,q1,alpha0,alpha1,alpha2):
		'''
			This correspond to the exact solution to the Grad-Shafranov equation in straight field line coordinates. 
			The equation reads:
			(PsiEdge/B0)**2 * R*(d(1/R*dPsiNorm/dR)/dR +d(dPsiNorm/dz)/dz)+1/2*(R**2 Beta'+R0**2*Y')
			
			where Beta = mu0/B0**2, Y=F**2/(R0*B0)**2, PsiNorm = Psi/PsiEdge and the prime is the derivative with respect to PsiNorm. Using the profiles
			
			Beta = 0.5*(alpha0 - alpha1 * PsiNorm)
			Y = alpha2- alpha3 * PsiNorm
			
			The solution reads
			
			Psi = (B0/(2*PsiEdge))**2 * (Z**2 *alpha3 *R0**2 + alpha1/4 * (R**2-R0**2)**2)
			
			This solution is useful because it can easily be parametrized in flux coordinates, and transformed into straight field line coordinates. The calculation is given in the Mathematica file Solovev.nb
			
			NOTE: I'm not 100% sure that the metric quantities are calculated correctly. I'm sure some of them are, particularly the ones used in the Grad-Shafranov equation (since it's fulfilled to machine precision).
		'''
		self.R0 = R0
		self.B0 = B0
		self.alpha0 = alpha0
		self.alpha1 = alpha1
		self.alpha2 = alpha2
		
		self.PsiEdge = B0*R0**2.*(alpha1*(q1**2.-q0**2.)/(16.*q1**2.-4.))**0.5
		self.alpha3 = alpha2/alpha1 * 4.*self.PsiEdge**2./(q0*B0*R0**2)**2.

		self.a = 4.*self.PsiEdge/(B0*R0**2.)*(1./alpha1)**0.5
		self.b = self.PsiEdge/(B0*R0)*(4./self.alpha3)**0.5
		

		self.gamma = 5./3. 
		self.mu0 = 4.*np.pi*1.0E-07


	def BuildInGrid(self,grid):
	
		#Variables
		#===================================		
		self.q = self.B0*self.R0*self.a*self.b*(self.alpha2-self.alpha3*grid.S**2.)**0.5/(4.*self.PsiEdge*(1.-self.a**2.*grid.S**2.)**0.5)
		self.F = (self.alpha2-self.alpha3*grid.S**2.)**0.5
		self.P = 0.5*(self.alpha0-self.alpha1*grid.S**2.)
		self.g = 2.*self.PsiEdge/(self.B0*self.R0**2.)
		self.h = grid.S
		self.rho = np.ones_like(grid.S)				#Constant density
		#self.rho = self.P/self.P[0]					#Constant temperature
		
		self.R = np.zeros((grid.Ntheta,grid.N+2))
		self.Z = np.zeros((grid.Ntheta,grid.N+2))
		
		self.dRds = np.zeros((grid.Ntheta,grid.N+2))
		self.dRdu = np.zeros((grid.Ntheta,grid.N+2))
		self.dZds = np.zeros((grid.Ntheta,grid.N+2))
		self.dZdu = np.zeros((grid.Ntheta,grid.N+2))
		
		self.dR2ds = np.zeros((grid.Ntheta,grid.N+2))
		self.dR2du = np.zeros((grid.Ntheta,grid.N+2))
		
		self.dB2ds = np.zeros((grid.Ntheta,grid.N+2))
		self.dB2du = np.zeros((grid.Ntheta,grid.N+2))
		
		self.dg11ds = np.zeros((grid.Ntheta,grid.N+2))
		self.dg11du = np.zeros((grid.Ntheta,grid.N+2))
		self.dg12ds = np.zeros((grid.Ntheta,grid.N+2))
		self.dg12du = np.zeros((grid.Ntheta,grid.N+2))
		self.dg22ds = np.zeros((grid.Ntheta,grid.N+2))
		self.dg22du = np.zeros((grid.Ntheta,grid.N+2))
		
		for i in range(grid.Ntheta):
		
			self.R[i] = ((1.-self.a**2.*grid.S**2.)/(1.+ self.a*grid.S*np.cos(grid.theta[i])))**0.5
			self.Z[i] = self.b*grid.S*np.sin(grid.theta[i])*(1.-self.a**2.*grid.S**2.)**0.5/(1.+self.a*grid.S*np.cos(grid.theta[i]))/self.R0
			
			self.dRds[i] = -(self.a*(2*grid.S*self.a + np.cos(grid.theta[i]) + grid.S**2*self.a**2*np.cos(grid.theta[i])))/(2.*np.sqrt((1 - grid.S**2*self.a**2)/(1 + grid.S*self.a*np.cos(grid.theta[i])))*(1 + grid.S*self.a*np.cos(grid.theta[i]))**2)
			self.dRdu[i]   = (grid.S*self.a*np.sqrt((1 - grid.S**2*self.a**2)/(1 + grid.S*self.a*np.cos(grid.theta[i])))*np.sin(grid.theta[i]))/(2 + 2*grid.S*self.a*np.cos(grid.theta[i]))
			self.dZds[i]   = -((self.b*(-1 + 2*grid.S**2*self.a**2 + grid.S**3*self.a**3*np.cos(grid.theta[i]))*np.sin(grid.theta[i]))/(np.sqrt(1 - grid.S**2*self.a**2)*self.R0*(1 + grid.S*self.a*np.cos(grid.theta[i]))**2))
			self.dZdu[i]   = (grid.S*np.sqrt(1 - grid.S**2*self.a**2)*self.b*(grid.S*self.a + np.cos(grid.theta[i])))/(self.R0*(1 + grid.S*self.a*np.cos(grid.theta[i]))**2)
			
			self.dR2ds[i]  = -((self.a*(2*grid.S*self.a + np.cos(grid.theta[i]) + grid.S**2*self.a**2*np.cos(grid.theta[i])))/(1 + grid.S*self.a*np.cos(grid.theta[i]))**2)
			self.dR2du[i]  = (grid.S*self.a*(1 - grid.S**2*self.a**2)*np.sin(grid.theta[i]))/(1 + grid.S*self.a*np.cos(grid.theta[i]))**2
		
			self.dg11ds[i] = (8*grid.S*self.a**4*self.R0**2 + self.a**3*self.R0**2*np.cos(grid.theta[i])*(4 + 20*grid.S**2*self.a**2 + 8*grid.S**4*self.a**4 + np.cos(grid.theta[i])*(2*grid.S*self.a*(1 + 16*grid.S**2*self.a**2 + 7*grid.S**4*self.a**4) + np.cos(grid.theta[i])*(-3 + grid.S**2*self.a**2 + 27*grid.S**4*self.a**4 + 7*grid.S**6*self.a**6 + grid.S*self.a*(-3 + 3*grid.S**2*self.a**2 + 7*grid.S**4*self.a**4 + grid.S**6*self.a**6)*np.cos(grid.theta[i])))) + 4*self.a*self.b**2*(-6*grid.S*self.a + 17*grid.S**3*self.a**3 - 6*grid.S**5*self.a**5 + (-4 + 8*grid.S**2*self.a**2 + 6*grid.S**4*self.a**4 - 3*grid.S**6*self.a**6)*np.cos(grid.theta[i]) + grid.S**3*self.a**3*(1 + 2*grid.S**2*self.a**2 + grid.S**3*self.a**3*np.cos(grid.theta[i]))*np.cos(2.*grid.theta[i]))*np.sin(grid.theta[i])**2)/(4.*(-1 + grid.S**2*self.a**2)**2*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**5)
			
			self.dg11du[i] = (-2*grid.S*self.a*(8*(7 - 27*grid.S**2*self.a**2 + 26*grid.S**4*self.a**4 + grid.S**6*self.a**6)*self.b**2 + self.a**2*(-17 + 38*grid.S**2*self.a**2 + 7*grid.S**4*self.a**4)*self.R0**2)*np.sin(grid.theta[i]) - 2*(8*(2 - 8*grid.S**2*self.a**2 + 3*grid.S**4*self.a**4 + 10*grid.S**6*self.a**6)*self.b**2 + self.a**2*(-4 + grid.S**2*self.a**2 + 30*grid.S**4*self.a**4 + grid.S**6*self.a**6)*self.R0**2)*np.sin(2.*grid.theta[i]) + grid.S*self.a*(-2*(8*(-1 + (grid.S*self.a + grid.S**3*self.a**3)**2)*self.b**2 + self.a**2*(-1 + 6*grid.S**2*self.a**2 + 7*grid.S**4*self.a**4)*self.R0**2)*np.sin(3.*grid.theta[i]) - grid.S*self.a**3*(8*grid.S**2*self.b**2 + (self.R0 + grid.S**2*self.a**2*self.R0)**2)*np.sin(4.*grid.theta[i])))/(32.*(-1 + grid.S**2*self.a**2)*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**5)
			
			self.dg22du[i] = -(grid.S**2*(-1 + grid.S**2*self.a**2)*((32*(-1 + 3*grid.S**2*self.a**2)*self.b**2 + self.a**2*(8 + 9*grid.S**2*self.a**2)*self.R0**2)*np.cos(grid.theta[i]) + grid.S*self.a*(16*(-1 + 4*grid.S**2*self.a**2)*self.b**2 + 14*self.a**2*self.R0**2 + 2*(8*self.b**2 + self.a**2*self.R0**2)*np.cos(2.*grid.theta[i]) - grid.S*self.a**3*self.R0**2*np.cos(3.*grid.theta[i])))*np.sin(grid.theta[i]))/(16.*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**5)
			
			self.dg22ds[i] = (grid.S*(32*(1 + grid.S**2*self.a**2 - 7*grid.S**4*self.a**4)*self.b**2 - self.a**2*(-8 + 17*grid.S**2*self.a**2 + grid.S**4*self.a**4)*self.R0**2 + 2*grid.S*self.a*(-8*(-9 + 20*grid.S**2*self.a**2 + 4*grid.S**4*self.a**4)*self.b**2 + self.a**2*(1 - 5*grid.S**2*self.a**2)*self.R0**2)*np.cos(grid.theta[i]) - 8*(4*(-1 + 3*grid.S**2*self.a**2 + grid.S**4*self.a**4)*self.b**2 + self.a**2*(1 - 2*grid.S**2*self.a**2)*self.R0**2)*np.cos(2.*grid.theta[i]) + grid.S*self.a*(-2*(8*self.b**2 + self.a**2*(1 - 5*grid.S**2*self.a**2)*self.R0**2)*np.cos(3.*grid.theta[i]) + grid.S*self.a**3*(1 + grid.S**2*self.a**2)*self.R0**2*np.cos(4.*grid.theta[i]))))/(32.*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**5)
			
			self.dg12ds[i] = -(((8*(-1 + 8*grid.S**2*self.a**2 + 5*grid.S**4*self.a**4)*self.b**2 + self.a**2*(2 + 7*grid.S**2*self.a**2)*self.R0**2)*np.cos(grid.theta[i]) + grid.S*self.a*(4*(-1 + 18*grid.S**2*self.a**2 + grid.S**4*self.a**4)*self.b**2 + self.a**2*(7 + grid.S**2*self.a**2)*self.R0**2 + (4*(3 + 2*grid.S**2*self.a**2 + grid.S**4*self.a**4)*self.b**2 + self.a**2*(-1 + grid.S**2*self.a**2)*self.R0**2)*np.cos(2.*grid.theta[i]) - grid.S*self.a**3*self.R0**2*np.cos(3.*grid.theta[i])))*np.sin(grid.theta[i]))/(8.*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**5)
		
			self.dg12du[i] = (grid.S*(-2*grid.S*self.a*(4*(-10 + 21*grid.S**2*self.a**2 + 6*grid.S**4*self.a**4)*self.b**2 + self.a**2*(15 + 19*grid.S**2*self.a**2)*self.R0**2)*np.cos(grid.theta[i]) + 4*(4*(2 - 7*grid.S**2*self.a**2 + 3*grid.S**4*self.a**4)*self.b**2 - self.a**2*(2 + grid.S**2*self.a**2 + grid.S**4*self.a**4)*self.R0**2)*np.cos(2.*grid.theta[i]) + grid.S*self.a*(20*grid.S*self.a*(4 - 9*grid.S**2*self.a**2)*self.b**2 - 5*grid.S*self.a**3*(9 + grid.S**2*self.a**2)*self.R0**2 + 2*(4*(-2 + grid.S**2*self.a**2 + 2*grid.S**4*self.a**4)*self.b**2 + self.a**2*(-1 + 3*grid.S**2*self.a**2)*self.R0**2)*np.cos(3.*grid.theta[i]) + grid.S*self.a**3*(self.R0**2 + grid.S**2*(4*self.b**2 + self.a**2*self.R0**2))*np.cos(4.*grid.theta[i]))))/(32.*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**5)
			
			self.dB2du[i]  = (grid.S*((4*self.a**2*self.PsiEdge**2)/self.b**2 + (self.a**4*(self.alpha2 - grid.S**2*self.alpha3)*self.B0**2*self.R0**2)/(-1 + grid.S**2*self.a**2) + (32*(-1 + grid.S**2*self.a**2)**2*self.PsiEdge**2)/(self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**3) + (4*(-1 + grid.S**2*self.a**2)*self.PsiEdge**2*(8*self.b**2 + self.a**2*self.R0**2))/(self.b**2*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**2))*np.sin(grid.theta[i]))/(self.a**3*self.B0**2)
			
			self.dB2ds[i]  = (2*grid.S*(self.a**2*self.alpha2 - self.alpha3)*self.R0**2)/(-1 + grid.S**2*self.a**2)**2 + (self.a*((-4*self.PsiEdge**2)/self.b**2 + (self.a**2*(self.alpha2 + grid.S**2*self.a**2*self.alpha2 + grid.S**2*(-3 + grid.S**2*self.a**2)*self.alpha3)*self.B0**2*self.R0**2)/(-1 + grid.S**2*self.a**2)**2)*np.cos(grid.theta[i]) + (4*self.PsiEdge**2*(2*grid.S*self.a + np.cos(grid.theta[i]) + grid.S**2*self.a**2*np.cos(grid.theta[i]))*(self.a*(8*grid.S**2*self.b**2 + self.R0**2) + grid.S*(8*self.b**2 + self.a**2*self.R0**2)*np.cos(grid.theta[i])))/(self.b**2*self.R0**2*(1 + grid.S*self.a*np.cos(grid.theta[i]))**3))/(self.a**2*self.B0**2)
			
			
		self.dqds   = (grid.S*(self.a**3*self.alpha2 - self.a*self.alpha3)*self.b*self.B0*self.R0)/(4.*(1 - grid.S**2*self.a**2)**1.5*np.sqrt(self.alpha2 - grid.S**2*self.alpha3)*self.PsiEdge)
		self.d2qds  = (self.a*(self.a**2*self.alpha2 - self.alpha3)*(self.alpha2 + 2*grid.S**2*self.a**2*self.alpha2 - 3*grid.S**4*self.a**2*self.alpha3)*self.b*self.B0*self.R0)/(4.*(1 - grid.S**2*self.a**2)**2.5*(self.alpha2 - grid.S**2*self.alpha3)**1.5*self.PsiEdge)
		self.dPds   = -self.alpha1*grid.S
		self.d2Pds  = -self.alpha1
		self.dFds   = -self.alpha3*grid.S/(self.alpha2-self.alpha3*grid.S**2.)**0.5
		
		self.R2     = self.R*self.R
		self.g11    = self.dRds**2 + self.dZds**2
		self.g12    = self.dRds*self.dRdu + self.dZds*self.dZdu
		self.g22    = self.dRdu**2 + self.dZdu**2
		
		self.B2     = self.F**2.*(self.q**2.*self.R2+self.g22)/(self.q**2.*self.R2**2)
		
		imag = 1.0j
		#===================================


	def plot(self,grid,show=True):

		ax2 = plt.subplot(231)
		plt.plot(grid.S,self.g*self.h*self.R0**2*self.B0)
		plt.tick_params('x', labelbottom=True)
		plt.xlabel('s')
		plt.ylabel(r'$f$ [Wb]')
		plt.grid()

		ax3 = plt.subplot(232)
		plt.plot(grid.S,self.F*self.R0*self.B0)
		plt.tick_params('x', labelbottom=True)
		plt.xlabel('s')
		plt.ylabel(r'$F$ [$Wb/m^2$]')
		plt.grid()
	
		ax4 = plt.subplot(234)
		plt.plot(grid.S,self.B0**2.*self.P/self.mu0/1000.)
		plt.tick_params('x', labelbottom=True)
		plt.xlabel('s')
		plt.ylabel('P [kPa]')
		plt.grid()
	
		ax5 = plt.subplot(235)
		plt.plot(grid.S,self.q)
		plt.ylabel(r'q')
		plt.xlabel('s')
		plt.grid()
		
		ax6 = plt.subplot(233)
		for i in range(0,grid.N+2,10):
			plt.plot(self.R[:,i]*self.R0,self.Z[:,i]*self.R0)
		if i != grid.N+1:
			plt.plot(self.R[:,grid.N+1]*self.R0,self.Z[:,grid.N+1]*self.R0, linewidth=2)
			
		plt.xlabel(r'R [m]')
		plt.ylabel(r'Z [m]')
		plt.grid()
		plt.gca().set_aspect('equal')

		
		
		##Check Grad Shafranov
		LHS = self.dFds*self.g22/(self.q*self.R2)+self.F*self.dg22ds/(self.q*self.R2)-self.F*self.g22*self.dqds/(self.q**2*self.R2)-self.F*self.g22*self.dR2ds/(self.q*self.R2**2.)-self.F/self.q*(self.dg12du/self.R2-self.g12*self.dR2du/self.R2**2.)
		RHS = -self.q*self.R2/self.F*(self.F*self.dFds/self.R2+self.dPds)
		
		x = 30

		plt.figure(2)
		plt.plot(grid.S,(RHS[x]-LHS[x])/max(abs(LHS[x])), 'o-', label=r'$\frac{|J^3-\Delta^*\psi_p/R^2|}{max(\Delta^*\psi_p/R^2)}$')
		plt.title('Check Grad-Shafranov equation')
		plt.grid()
		

		plt.legend()
		if show: plt.show()



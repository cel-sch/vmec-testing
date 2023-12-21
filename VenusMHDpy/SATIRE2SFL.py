import h5py
from scipy.io import netcdf
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from .library import Fsincos,Fsin,Fcos,SplineDer1,SplineDer2


#***************************************************************************************************************#
#                                                                                                               #
#                                                SATIRE2SFL                                                     #
#                                                                                                               #
#***************************************************************************************************************#


class SATIRE2SFL:
	def __init__(self,woutfile):
		'''
			This file maps a VMEC equilibria into a Straigth Field Line right-handed coordinate system. The result is 
			written into an .h5 file 'eq.sfl.h5', which is later read by the stability solver. 

			The radial variable is defined to be an evenly space grid s2 = [0.,1.] interval. Here, an evenly spaced grid
			represents the normalized poloidal (LRFP=True) OR torodial (LRFP=False) flow, whatever was chosen in VMEC. 
			The stability code will later use the same variable as the radial variable by default. 

			The derivatives are taken with respect to rho = s = sqrt(s2).

			Author: Guillermo Bustos Ramirez
			Last modified: 08/April/2022 
		'''


		#Read the variables from the VMEC wout_.nc file
		#=====================================================================================================
		print ('--------------------------------SATIRE2SFL----------------------------------')
		print ('Reading SATIRE (VMEC-FLOW) equilibrium...')
		time0 = time.time()

		Wout = netcdf.netcdf_file(woutfile,'r',mmap=False)
		LASYM = Wout.variables['lasym__logical__'][()].copy()	#Equlibrium violates stellarator symmetry. In tokamaks, this simplifies as up-down symmetry.
		LRFP  = Wout.variables['lrfp__logical__'][()].copy()	#Wheather equilibrium uses normalised toroidal (F) or poloidal (T) flux as a radial variable.
		xm = Wout.variables['xm'][()].copy()					#Poloidal modes allowed
		xn = Wout.variables['xn'][()].copy()					#Toroidal modes allowed
		rmnc = Wout.variables['rmnc'][()].copy()				#Cosine component of cylindrical variable R
		zmns = Wout.variables['zmns'][()].copy()				#Sine component of cylindrical variable Z
		lmns = Wout.variables['lmns'][()].copy()				#Sine component of VMEC's lambda stream function (half grid)

		if LASYM: 
			print ('Equilibrium calculated WITHOUT assuming stellarator symmetry')
			
			rmns = Wout.variables['rmns'][()].copy()	#Sine component of cylindrical variable R.
			zmnc = Wout.variables['zmnc'][()].copy()	#Cosine component of cylindrical variable Z.
			lmnc = Wout.variables['lmnc'][()].copy()	#Cosine component of VMEC's lambda stream function (half grid)
		else:
			print ('Equilibrium calculated assuming stellarator symmetry')


		ns = Wout.variables['ns'][()].copy()					#Radial discretisation: Number of closed flux surfaces
		PsiP_VMEC = Wout.variables['chi'][()].copy()			#VMEC's poloidal flux. PsiP = PsiP_VMEC/(2.*pi)
		dPsiP_VMEC = Wout.variables['chipf'][()].copy()			#VMEC's derivative of poloidal flux. dPsiP = dPsiP_VMEC/(2.*pi)
		PsiT_VMEC = Wout.variables['phi'][()].copy()			#VMEC's toroidal flux. PsiT = PsiT_VMEC/(2.*pi)
		dPsiT_VMEC = Wout.variables['phipf'][()].copy()			#VMEC's derivative of toroidal flux. dPsiT = dPsiT_VMEC/(2.*pi)


		q = Wout.variables['q_factor'][()].copy()			#Safety factor profile
		Pbar = Wout.variables['presf'][()].copy()			#Pressure profile (Pbar(s)). Pbar does NOT include the effects of rotation, so the pressure profile should be P(s,u) = Pbar(s)*exp(U(s)*(R(s,u)**2-R0**2))
		F  = Wout.variables['bvco'][()].copy()				#Poloidal average of toroidal current (half grid). Function F in the magnetic field B = GradPsi X Grad Phi + F * GradPhi
		AM = Wout.variables['am'][()].copy()				#Polynomial coefficients of pressure profile (P_(s)). In VMEC, P(s,u) = P_(s)*exp(U(s)*R(s,u)), therefore, Pbar(s) = P_(s)*exp*(U(s)*R0).
		AT = Wout.variables['at'][()].copy()				#Polynomial coefficients of normalised temperature profile (T(s))
		AH = Wout.variables['ah'][()].copy()				#Polynomial coefficients of normalised rotation profile (Omega(s)).
		self.B0  = Wout.variables['b0'][()].copy()			#Magnetic field on axis
		self.beta0 = Wout.variables['betatotal'][()].copy()	#Plasma beta on axis
		
		try:
			self.M02 = Wout.variables['machsq'][()].copy()			#Mach number on axis squared. M02 = M0^2 = mi*Omega0^2*R0^2/(2*T0)
			FlowVersion = True
		except KeyError:
			self.M02 = 0.
			FlowVersion = False

		#Close the file
		Wout.close()

		print ('Equilibrium read in %.4f s'%(time.time()-time0))
		#=====================================================================================================




		#===============================GRID=====================================#
		#Straight Field Line poloidal angle (equidistant)
		self.ntheta = 300
		theta = np.linspace(0.,2.*np.pi,self.ntheta,endpoint=False) 

		#Poloidal flux
		PsiP = -PsiP_VMEC/(2.*np.pi)

		#Toroidal flux
		PsiT = -PsiT_VMEC/(2.*np.pi)

		#VMEC normalised radial variable (Can be poloidal or toroidal flux)
		s2 = np.linspace(0.,1.,ns)
		self.s  = np.sqrt(s2)
		self.h  = self.s
	
		#Normalised toroidal/poloidal flux in half grid
		s2_half = s2[:-1]+np.diff(s2)*0.5
		#========================================================================#





		time0 = time.time()
		#============================GEOMETRY====================================#
		#Can be made more efficient by getting rid of the loops and performing vector opeartions.
		#It's fast enough, and has to be performed only once, so it doesn't matter too much, 

		print ('Mapping surfaces to Straight Field Line coordinates...')
		self.R = np.zeros((self.ntheta,ns))
		self.Z = np.zeros((self.ntheta,ns))
		LAMBDA = np.zeros((self.ntheta,ns-1))

		#The stream function LAMBDA gives the transformation to SFL coordinates. 
		#Given in half grid. Need interpolation (and extrapolation in the two boundary points)
		#to obtain in full grid.
		#---------------------------------------------------------------------------------------
		for i in range(ns-1):
			if LASYM:
				LAMBDA[:,i] = Fsincos(theta,0.,xm,xn,lmnc[i+1],lmns[i+1],False)
			else:
				LAMBDA[:,i] = Fsin(theta,0.,xm,xn,lmns[i+1],False)
			
		dummy = interpolate.interp1d(s2_half,LAMBDA,kind='cubic',axis=1,fill_value='extrapolate')
		LAMBDAf = dummy(s2)   
		#---------------------------------------------------------------------------------------

		for i in range(ns):
			#Reconstruction of the coordinates from Fourier coefficients
			if LASYM:
				R_ = Fsincos(theta,0.,xm,xn,rmnc[i],rmns[i],False)
				Z_ = Fsincos(theta,0.,xm,xn,zmnc[i],zmns[i],False)
			else:
				R_ = Fcos(theta,0.,xm,xn,rmnc[i],False)
				Z_ = Fsin(theta,0.,xm,xn,zmns[i],False)

			#Interpolation into equidistant straight field line poloidal angle.
			dummy = interpolate.splrep(theta+LAMBDAf[:,i],R_)
			self.R[:,i] = interpolate.splev(theta,dummy)
			dummy = interpolate.splrep(theta+LAMBDAf[:,i],Z_)
			self.Z[:,i] = interpolate.splev(theta,dummy)

		#========================================================================#


			

		#============================PROFILES====================================#
		print ('Extracting profiles...')
		#Physical constants
		#-------------------------------------------------
		#Major radius
		self.R0 = rmnc[0][0]
		self.gamma = 5./3. 
		self.mu0 = 4.*np.pi*1.0E-07
		self.kappa = 0.5*np.sqrt(np.pi)
		self.Kpar = 0.
		
		self.c0 = 1.
		self.c1 = 1.
		self.c2 = 1.
		self.c3 = 1.
		
		
		#Toroidal magnetic field function F
		#-------------------------------------------------
		dummy = interpolate.splrep(s2_half,F[1:])
		self.F = -interpolate.splev(s2,dummy)


		#Safety factor profile. Note the change in the sign
		#due to VMEC using left-hand coordinate system.
		#-------------------------------------------------
		self.q = -q


		#Normalised temperature profile (T/T0). But note that
		#the absolute value of the array will be the same as
		#the one given as input, BUT! VMEC always uses the 
		#normalised temperature as T --> T/T[0]. Therefore,
		#T is normalised with respect to its value on axis.
		#-------------------------------------------------
		if FlowVersion:
			self.T = np.polyval(AT[::-1],s2)
			self.T = self.T/self.T[0]
			if self.T[-1]<0: self.T[-1]= -self.T[-1]


		#Normalised rotation profile Omega. As with the 
		#temperature, VMEC uses the normalised rotation
		#profile as Omega --> Omega/Omega[0]. So we need
		#to normalise the extracted profile.
		#-------------------------------------------------
		if FlowVersion:
			self.Omega = np.polyval(AH[::-1],s2)
			self.Omega = self.Omega/self.Omega[0]
		else:
			self.Omega = np.zeros_like(s2)
		
		
		#VMEC U profile. Equilibrium in VMEC is perturbed by rotation by
		#modifying the pressure profile as 
		#P(s,u) = P_(s)*exp(U(s)*(R(s,u)**2-R0**2.), where
		#U(s) = mi*Omega**2./(4*T) = M02*(Omega/Omega0)**2./(2.*R0**2*(T/T0)).
		#Since we only know the normalised Omega and T, then
		#the second expression is used, where 
		#M02=mi*Omega0**2*R0**2/(2*T0) is the Mach number squared
		#-------------------------------------------------
		if FlowVersion:
			self.U = self.M02*self.Omega**2./(2.*self.R0**2.*self.T)
			#self.U[-1] = 128.*self.M02/(np.pi*self.R0)**2.
		else:
			self.U = np.zeros_like(s2)
		
		#Pressure profile P_(s) = P(s)*exp(-U(s)*R0**2), where
		#P(s) is the pressure profile in the absence of rotation.
		#Actual pressure profile is 
		#Prot(s,u) = P(s)*exp(U(s)*(R(s,u)**2-R0**2)) = P_(s)*exp(U(s)*R(s,u)**2)
		#-------------------------------------------------
		P_   = np.polyval(AM[::-1],s2)
		self.P    = P_*np.exp(self.U*self.R0**2.)
		self.Prot = P_*np.exp(self.U*self.R**2.)
		self.P0   = self.P[0]
			
		#Density is rho(s) = P(s)/T(s) in the absence of rotation.
		#Note that with the info provided by VMEC, we can't 
		#know the absolute value of the density, so only the
		#normalised value is calculated.
		#rho(s) = P(s)/2T = P(s)/(2T0*T/T0) = P(s)*rho0/(T*P0/T0) --> rho(s)/rho0 = P(s)/((T/T0)*P0)
		#-------------------------------------------------
		if FlowVersion:
			self.rho    = self.P/(self.T*self.P0)
			self.rhorot = self.Prot/(self.T*self.P0)
		else:
			#Constant density
			#self.rho = np.ones_like(s2)		
			#self.T = self.P/(self.P0*self.rho)
			
			#Constant temperature
			self.T = np.ones_like(s2)
			self.rho = self.P/(self.T*self.P0)

			self.rhorot = self.Prot/(self.T*self.P0)


		#Thermal ion velocity.
		#Uthi = (Ti)**0.5. But since we don't know T0, we need
		#to provide the normalised value Uthi/vA, where
		#	Va = B0/(mu0*rho0)**0.5 is the Alfven velocity on axis
		#We then have
		#Uthi/vA = ((T/T0)*mu0*P0/(2.*B0**2.))**1/2
		#-------------------------------------------------
		self.Uthi = (self.T*self.mu0*self.P0/(2.*self.B0**2.))**0.5 


		#Magnetic flux derivatives
		#-------------------------------------------------
		if LRFP:
			PsiP_prime = 2.*self.s*PsiP[-1]
			PsiT_prime = self.q*2.*self.s*PsiP[-1]
			self.g = 2.*PsiP[-1]*np.ones_like(self.s)
		else: 
			PsiP_prime = 2.*self.s*PsiT[-1]/self.q
			PsiT_prime = 2.*self.s*PsiT[-1]
			self.g = 2.*PsiT[-1]/self.q
		#========================================================================#
		print ('Calculations completed in %.4f s'%(time.time()-time0))
		print ('-----------------------------------------------------------------------------\n')



	def Writeh5(self,eqFile):
		#========================WRITE h5 FILE===================================#
		time0 = time.time()
		print ('Writing equilibrium into h5 file...')

		h5 = h5py.File(eqFile,'w')

		#Geometry
		geometry = h5.create_group('geometry')

		geometry.create_dataset('R',data=self.R)
		geometry.create_dataset('Z',data=self.Z)

		#Normalisation
		normalisation  = h5.create_group('normalisation')

		normalisation.create_dataset('R0',data=self.R0)
		normalisation.create_dataset('B0',data=self.B0)
		normalisation.create_dataset('P0',data=self.P0)
		normalisation.create_dataset('M02',data=self.M02)

		#Profiles
		profiles = h5.create_group('profiles')

		profiles.create_dataset('F',data=self.F)
		profiles.create_dataset('g', data=self.g)
		profiles.create_dataset('h', data=self.s)
		profiles.create_dataset('q',data=self.q)
		profiles.create_dataset('P',data=self.P)
		profiles.create_dataset('Prot',data=self.Prot)
		profiles.create_dataset('rho',data=self.rho)
		profiles.create_dataset('rhorot',data=self.rhorot)
		profiles.create_dataset('T',data=self.T)
		profiles.create_dataset('U',data=self.U)
		profiles.create_dataset('s',data=self.s)

		h5.close()
		print ('h5 file written in %.4f s'%(time.time()-time0))
		print ('-------------------------------------------------------------')
		#========================================================================#
		

	def Normalise(self):
		
		
		#Mass density rho is already normalised with respect to its value on axis rho0.
		#Ion thermal velocity Uthi is already normalised with respect to Alfven velocity on axis.
		self.Omega = self.Omega*(self.M02*self.mu0*self.P0)**0.5/self.B0
		self.Prot  = self.mu0*self.Prot/self.B0**2.
		self.F = self.F/(self.R0*self.B0)
		self.g = self.g/(self.R0**2.*self.B0)
		self.P = self.mu0*self.P/self.B0**2.
		self.T = self.mu0*self.P0*self.T/(2.*self.B0**2.)
		self.U = self.U*self.R0**2.
		self.R = self.R/self.R0
		self.Z = self.Z/self.R0
		

	def ChangeGrid(self,S_new):
	
		dummy  = interpolate.splrep(self.s,self.q)
		self.q = interpolate.splev(S_new,dummy)

		dummy  = interpolate.splrep(self.s,self.F)
		self.F = interpolate.splev(S_new,dummy)
		
		dummy  = interpolate.splrep(self.s,self.g)
		self.g = interpolate.splev(S_new,dummy)
		
		dummy  = interpolate.splrep(self.s,self.h)
		self.h = interpolate.splev(S_new,dummy)

		dummy  = interpolate.splrep(self.s,self.P)
		self.P = interpolate.splev(S_new,dummy)
		
		dummy  = interpolate.splrep(self.s,self.rho)
		self.rho = interpolate.splev(S_new,dummy)
		
		dummy  = interpolate.splrep(self.s,self.T)
		self.T = interpolate.splev(S_new,dummy)

		dummy  = interpolate.splrep(self.s,self.U)
		self.U = interpolate.splev(S_new,dummy)
		
		dummy  = interpolate.splrep(self.s,self.Omega)
		self.Omega = interpolate.splev(S_new,dummy)

		dummy  = interpolate.splrep(self.s,self.Uthi)
		self.Uthi = interpolate.splev(S_new,dummy)
		
		dummy  = interpolate.interp1d(self.s,self.Prot,kind='cubic',axis=1)
		self.Prot = dummy(S_new)  

		dummy  = interpolate.interp1d(self.s,self.rhorot,kind='cubic',axis=1)
		self.rhorot = dummy(S_new)  

		dummy  = interpolate.interp1d(self.s,self.R,kind='cubic',axis=1)
		self.R = dummy(S_new)   

		dummy  = interpolate.interp1d(self.s,self.Z,kind='cubic',axis=1)
		self.Z = dummy(S_new)
		
		self.s = S_new

	def BuildInGrid(self,grid):
	
		#Variables
		#===================================
		#Profiles
		self.drhods        = SplineDer1(self.rho,grid,der_s=1)
		self.dTds          = SplineDer1(self.T,grid,der_s=1)
		self.dUds          = SplineDer1(self.U,grid,der_s=1)
		self.dOmegads      = SplineDer1(self.Omega,grid,der_s=1)
		self.dgds          = SplineDer1(self.g,grid,der_s=1)
		self.dhds          = SplineDer1(self.h,grid,der_s=1)
		self.dqds          = SplineDer1(self.q,grid,der_s=1)
		self.d2qds         = SplineDer1(self.q,grid,der_s=2)
		self.dPds          = SplineDer1(self.P,grid,der_s=1)
		self.d2Pds         = SplineDer1(self.P,grid,der_s=2) 
		self.dFds          = SplineDer1(self.F,grid,der_s=1)
		self.drhorotds     = SplineDer2(self.rhorot,grid,der_s=1,der_theta=0)
		self.drhorotdu     = SplineDer2(self.rhorot,grid,der_s=0,der_theta=1)
		self.d2rhorotds    = SplineDer2(self.rhorot,grid,der_s=2,der_theta=0)
		self.d2rhorotdu    = SplineDer2(self.rhorot,grid,der_s=0,der_theta=2)
		self.d2rhorotdsdu  = SplineDer2(self.rhorot,grid,der_s=1,der_theta=1)
		self.dProtds       = SplineDer2(self.Prot,grid,der_s=1,der_theta=0)
		self.dProtdu       = SplineDer2(self.Prot,grid,der_s=0,der_theta=1)
		self.d2Protds      = SplineDer2(self.Prot,grid,der_s=2,der_theta=0)
		self.d2Protdu      = SplineDer2(self.Prot,grid,der_s=0,der_theta=2)
		self.d2Protdsdu    = SplineDer2(self.Prot,grid,der_s=1,der_theta=1)
		
		#Geometry
		self.dRds          = SplineDer2(self.R,grid,der_s=1,der_theta=0)
		self.dRdu          = SplineDer2(self.R,grid,der_s=0,der_theta=1)
		self.dZds          = SplineDer2(self.Z,grid,der_s=1,der_theta=0)
		self.dZdu          = SplineDer2(self.Z,grid,der_s=0,der_theta=1)
		self.R2            = self.R*self.R
		self.dR2ds         = SplineDer2(self.R2,grid,der_s=1,der_theta=0)
		self.d2R2ds        = SplineDer2(self.R2,grid,der_s=2,der_theta=0)
		self.dR2du         = SplineDer2(self.R2,grid,der_s=0,der_theta=1)
		self.d2R2du        = SplineDer2(self.R2,grid,der_s=0,der_theta=2)
		self.d2R2dsdu      = SplineDer2(self.R2,grid,der_s=1,der_theta=1)
		self.g11           = self.dRds**2 + self.dZds**2
		self.dg11ds        = SplineDer2(self.g11,grid,der_s=1,der_theta=0)
		self.dg11du        = SplineDer2(self.g11,grid,der_s=0,der_theta=1)
		self.g12           = self.dRds*self.dRdu + self.dZds*self.dZdu
		self.dg12ds        = SplineDer2(self.g12,grid,der_s=1,der_theta=0)
		self.dg12du        = SplineDer2(self.g12,grid,der_s=0,der_theta=1)
		self.g22           = self.dRdu**2 + self.dZdu**2
		self.dg22ds        = SplineDer2(self.g22,grid,der_s=1,der_theta=0)
		self.dg22du        = SplineDer2(self.g22,grid,der_s=0,der_theta=1)
		
		#Fields
		self.B2            = self.F**2.*(self.q**2.*self.R2+self.g22)/(self.q**2.*self.R2**2)
		self.dB2ds         = SplineDer2(self.B2,grid,der_s=1,der_theta=0)
		self.dB2du         = SplineDer2(self.B2,grid,der_s=0,der_theta=1)
		self.J2 = -((self.dFds*self.F)/(self.g*self.h*self.q*self.R2))
		#self.J3 = self.F/(self.q*self.R2*self.g*self.h)*(self.dFds*self.g22/(self.q*self.R2)+self.F*self.dg22ds/(self.q*self.R2)-self.F*self.g22*self.dqds/(self.q**2*self.R2)-self.F*self.g22*self.dR2ds/(self.q*self.R2**2.)-self.F/self.q*(self.dg12du/self.R2-self.g12*self.dR2du/self.R2**2.))
		self.J3 =  self.F/(self.q*self.R2*self.g*self.h)*( -self.q*self.dFds-(self.q*self.R2/self.F)*(self.dProtds-self.Prot*self.U*self.dR2ds))
		
		
		imag = 1.0j
		#===================================


	def plot(self,grid,show=True):

		#ax1 = plt.subplot(251)
		#plt.plot(grid.S,self.F*self.R0*self.B0)
		#plt.tick_params('x', labelbottom=False)
		#plt.ylabel(r'F $[Wb/m]$')
		#plt.xlabel('s')
		#plt.grid()

		#ax2 = plt.subplot(252)
		#plt.plot(grid.S,self.g*self.h*self.R0**2.*self.B0)
		#plt.tick_params('x', labelbottom=False)
		#plt.xlabel('s')
		#plt.ylabel('f [Wb]')
		#plt.grid()
		
		ax3 = plt.subplot(241)
		plt.plot(grid.S,2.*self.B0**2.*self.T/(self.P0*self.mu0))
		plt.ylabel(r'$T/T_0$')
		plt.xlabel('s')
		plt.grid()

		ax4 = plt.subplot(242)
		plt.plot(grid.S,self.q)
		plt.ylabel(r'q')
		plt.xlabel('s')
		plt.grid()

		ax5 = plt.subplot(243)
		plt.plot(grid.S,self.B0**2.*self.P/self.mu0, label=r'$P(s)$')
		plt.plot(grid.S,self.B0**2.*self.Prot[0]/self.mu0, label=r'$P(s)e^{U(R^2-R_0^2)}$')
		plt.legend()
		plt.tick_params('x', labelbottom=False)
		plt.xlabel('s')
		plt.ylabel(r'$\bar{P}$ [Pa]')
		plt.grid()
	
		ax6 = plt.subplot(244)
		plt.plot(self.R[:,-1]*self.R0,self.Z[:,-1]*self.R0)
		plt.plot(self.R[0,0]*self.R0,self.Z[0,0]*self.R0,'+')
		plt.xlabel(r'R [m]')
		plt.ylabel(r'Z [m]')
		plt.gca().set_aspect('equal')
		plt.grid()
		
		ax7 = plt.subplot(245)
		plt.plot(grid.S,self.U/self.R0**2.)
		plt.ylabel(r'$U$ $[m^{-2}]$')
		plt.xlabel('s')
		plt.grid()
		
		ax8 = plt.subplot(246)
		plt.plot(grid.S,self.rho)
		plt.ylabel(r'$\bar{\rho}/\rho_0$')
		plt.xlabel('s')
		plt.grid()

		ax9 = plt.subplot(247)
		if self.M02 == 0:
			plt.plot(grid.S,np.zeros_like(grid.S))
		else:
			plt.plot(grid.S,self.Omega*self.B0/(self.M02*self.mu0*self.P0)**0.5)
		plt.ylabel(r'$\Omega/\Omega_0$')
		plt.xlabel('s')
		plt.grid()
		
		
		
		#Check Grad Shafranov
		#LHS = (SplineDer2(self.F*self.g22/(self.q*self.R2),grid,der_s=1,der_theta=0)-SplineDer2(self.F*self.g12/(self.q*self.R2),grid,der_s=0,der_theta=1))
		LHS = self.dFds*self.g22/(self.q*self.R2)+self.F*self.dg22ds/(self.q*self.R2)-self.F*self.g22*self.dqds/(self.q**2*self.R2)-self.F*self.g22*self.dR2ds/(self.q*self.R2**2.)-self.F/self.q*(self.dg12du/self.R2-self.g12*self.dR2du/self.R2**2.)
		RHS = -self.q*self.dFds-(self.q*self.R2/self.F)*(self.dPds+(self.R**2.-1.)*self.P*self.dUds)*np.exp(self.U*(self.R**2.-1.))
		#RHS = -self.q*self.dFds-(self.q*self.R2/self.F)*(self.dProtds-self.Prot*self.U*self.dR2ds)
		
		x = 0
		plt.figure(2)
		#plt.plot(grid.S,LHS[x], 'x-', label=r'$\Delta^*\psi_p (\mathcal{J}/R^2)$')
		#plt.plot(grid.S,RHS[x], 'o-', label=r'$J^\phi\mathcal{J}$')
		plt.plot(grid.S,(RHS[x]-LHS[x])/max(abs(LHS[x])), 'o-', label=r'$\frac{|J^\phi\mathcal{J}-\Delta^*\psi_p(\mathcal{J}/R^2)|}{max[\Delta^*\psi_p(\mathcal{J}/R^2)]}$')
		plt.grid()

		plt.legend()
		if show: plt.show()


if __name__ == '__main__':

	#Open the file
	try:
		woutfile = sys.argv[1]

		try:
			eqFile = sys.argv[2]
		except IndexError:
			eqFile = 'eq.SFL.h5'

		dummy = VMEC2SFL(woutfile)
		dummy.Writeh5(eqFile)

	except IndexError:
		None







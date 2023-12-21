from datetime import datetime
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigs
from .Bsplines import Bspline


class SolveScipy:
	def __init__(self,A,B,EV_guess=None,N_EV=1,maxiter=1000, which='LR', EValuesFile=None, EValuesID=0., EVectorsFile=None, grid=None, Model=None, ntor=None):
		'''
			Given the discretised general eigenvalue problem Ax = lambda*Bx, it finds the solution using the scypy.sparse.linalg.eigs function of the scipy library.
			
			Input:
					-A, B:          Matrices A and B of the general eigenvalue problem
					-EV_guess:      Eigenvalue EV_guess
					-N_EV:          Number of eigenvalues to obtain
					-maxiter:       Maximum number of iterations
					-Values_txt:    Name of the txt file where eigenvalues are to be stored. It appends the result to the existing file so that it's easier to recover eigenvalues for parameter scans
					-ValuesID:      Label for the simulation that is to be written in the .txt file
					-Vectors_h5:    Name of the .h5 file where eigenvectos are to be stored.
					-grid:          Grid object used to discretise the problem. 
					-Model:         MHD model used. Only needed for printing information
					-ntor:			Toroidal mode number. Only needed for printing information
					
		'''

		EVkind = 'Not applicable'
		if Model in ['IdealMHD','IdealMHDBussac','IdealMHDperp','IdealMHDcyl']:
			EVkind = '(gamma/Omega_A)^2'
		if Model in ['IdealMHDEuler','IdealMHDFlow-Euler','CASTOR']:
			EVkind = '(gamma/Omega_A)'

		#Solve
		#------------------------------------------------------------------------------------
		self.vals, self.vecs = eigs(A,k=N_EV,M=B,sigma=EV_guess,maxiter=maxiter, which=which)
		#self.vals, self.vecs = eigs(A,k=N_EV,M=B,sigma=EV_guess,maxiter=maxiter, which=which, v0=np.zeros(grid.N+1+max(grid.nu)))
		#------------------------------------------------------------------------------------
		
		#Write the .txt file with eigenvalues 
		self.grid = grid
		if EValuesFile != None:
			F = open(EValuesFile, 'a+')
			F.write('#--------------------------------- VENUS-MHDpy ---------------------------------# \n \n')
			F.write('Discretised using the %s model \n'%(Model))
			F.write('Normalised growth rates at the magnetic axis: %s \n'%(EVkind))
			F.write('Run on %s \n'%(datetime.now()))
			F.write('Nradial   = %i \n'%(self.grid.N))
			F.write('Npoloidal = %i \n'%(self.grid.Ntheta))
			F.write('n = %i \n'%(ntor))
			F.write('m = %i,%i \n'%(self.grid.Mmin,self.grid.Mmax))
			
			if EV_guess.imag == 0.:
				jsign = 1.
			else:
				jsign = EV_guess.imag/abs(EV_guess.imag)
				
			F.write('EV guess = %.6E + (%i)* %.6Ej \n'%(EV_guess.real,jsign,abs(EV_guess.imag)))
			F.write('Simulation ID: %s \n'%(str(EValuesID)))
			F.write('\n')
			F.write('Eigenvalue list: \n')
			for val in self.vals:
				jsign = val.imag/abs(val.imag)
				F.write(' %.6E + (%i)* %.6Ej \n'%(val.real,jsign,abs(val.imag)))
			F.write('#-------------------------------------------------------------------------------# \n')
			F.close()
		
		#Write the .h5 file with eigenvectors
		if EVectorsFile!=None:
			if grid==None:
				print ('Grid missing, cannot generate eigen vectors .h5 file...')
			else:
				h5 = h5py.File(EVectorsFile, 'w')
				
				#Eigenvalues
				EigenValues = h5.create_dataset('EigenValues',data=self.vals)
				
				#Grid data
				Grid = h5.create_group('grid')
				
				Grid.create_dataset('N',data=grid.N)
				Grid.create_dataset('nu',data=grid.nu)
				Grid.create_dataset('S',data=grid.S)
				Grid.create_dataset('Mtot',data=grid.Mtot)
				Grid.create_dataset('m',data=grid.m)
				Grid.create_dataset('knots',data=grid.knots)
				Grid.create_dataset('sk',data=grid.sk)
				
				#Perturbed quantities
				Variables = h5.create_group('variables')
				for i,nu in enumerate(grid.nu):
					dummy = np.zeros(shape=(grid.Mtot,grid.N+1+nu), dtype=complex)
					mixB = max(grid.nu)-nu
					  
					for j in range(grid.Mtot):
						dummy[j] = self.vecs[mixB+i*grid.dimension*grid.Mtot+j*grid.dimension:i*grid.dimension*grid.Mtot+(j+1)*grid.dimension,0]
						
					Variables.create_dataset('var'+str(i),data=dummy)


				h5.close()
	

	def PlotEigenValues(self, show=False):

		#Plot Eigenvalues
		#-------------------------------------------
		plt.figure()
		plt.plot(self.vals.real,self.vals.imag,'x')
		#plt.ylim((-1.E-03,1.E-03))
		plt.xlabel(r'$\mathcal{R}$')
		plt.ylabel(r'$\mathcal{I}$')
		plt.grid()
		if show: plt.show()
			
	def PlotEigenVectors(self, eq=None, PlotDerivatives=False, N_plot=0, show=True):
		'''
			Plots eigenvectors and first derivatives along with the safety factor.
			
			Input:
				-grid:              Grid object used to discretise the problem. 
				-eq:                Equilibrium object used to extract the safety factor.
				-PlotDerivatives:   Whether to additionally plotting the eigenfunctions, plotting the derivatives.
				-N_plot:            Eigenvectors corresponding to which eigenvalue to plot
				-show:              Show the figure
		'''
		
		#delta = 1.E+03
		#r = np.linspace(-delta,delta,10000)
		r = np.linspace(0.,1.,10000)
		#r = self.grid.S
		BspCalc = []

		for i,nu in enumerate(self.grid.nu):			
	
			#Create the Bspline arrays
			#-------------------------------------------------------------------------------------------------
			if nu not in BspCalc:
				BspCalc.append(nu)	
				vars(self)['Bsp'+str(nu)] = np.zeros(shape=(r.size, self.grid.N+1+nu))
				if PlotDerivatives: vars(self)['dBsp'+str(nu)] = np.zeros(shape=(r.size, self.grid.N+1+nu))
				
				for j in range(self.grid.N+1+nu):
					l = np.ones(len(r), dtype=int)*j-nu
					vars(self)['Bsp'+str(nu)][:,j] = Bspline(r,l,nu,self.grid,der=0)
					if PlotDerivatives: vars(self)['dBsp'+str(nu)][:,j] = Bspline(r,l,nu,self.grid,der=1)
			#-------------------------------------------------------------------------------------------------
		

			if PlotDerivatives:
				fig, vars()['ax'+str(i)] = plt.subplots(2, sharex=True)
				dlns = []
			else:
				fig, vars()['ax'+str(i)] = plt.subplots(1)
			 
			mixB = max(self.grid.nu)-nu
			lns = []
			for j in range(self.grid.Mtot):
				
				vec = self.vecs[mixB+i*self.grid.dimension*self.grid.Mtot+j*self.grid.dimension:i*self.grid.dimension*self.grid.Mtot+(j+1)*self.grid.dimension,N_plot].real
				mode = np.sum(vars(self)['Bsp'+str(nu)]*vec, axis=1)
	
				if PlotDerivatives: 
					dmode = np.sum(vars(self)['dBsp'+str(nu)]*vec, axis=1)
					lns_  = vars()['ax'+str(i)][0].plot(r,mode,'-', lw=2, label='m=%i'%(self.grid.m[j]))
					dlns_ = vars()['ax'+str(i)][1].plot(r,dmode,'-', lw=2, label='m=%i'%(self.grid.m[j]))
					lns += lns_
					dlns += dlns_
				else:
					lns_  = vars()['ax'+str(i)].plot(r,mode,'-', lw=2, label='m=%i'%(self.grid.m[j]))
					lns += lns_
									
		
			if PlotDerivatives:
                
	
				labs = [l.get_label() for l in lns]
				dlabs = [l.get_label() for l in dlns]

				vars()['ax'+str(i)][0].legend(lns, labs, loc='best')
				vars()['ax'+str(i)][1].legend(dlns, dlabs, loc='best')
				
				vars()['ax'+str(i)][0].grid()
				vars()['ax'+str(i)][1].grid()
				vars()['ax'+str(i)][0].tick_params(axis='y', labelcolor='red')
				vars()['ax'+str(i)][1].tick_params(axis='y', labelcolor='red')
				vars()['ax'+str(i)][1].set_xlabel(r'$s=\sqrt{\Psi/\Psi[1]}$')
				vars()['ax'+str(i)][0].set_ylabel(r'au', color='red')
				vars()['ax'+str(i)][1].set_ylabel(r'au', color='red')
				fig.suptitle('Variable %i'%(i+1))
				
				if eq != None:
					vars()['ax0'+str(i)] = vars()['ax'+str(i)][0].twinx()
					vars()['ax1'+str(i)] = vars()['ax'+str(i)][1].twinx()

					vars()['ax0'+str(i)].plot(self.grid.S,eq.q, '-', color='blue', lw=2, label='q')
					vars()['ax1'+str(i)].plot(self.grid.S,eq.q, '-', color='blue', lw=2, label='q')
					
					vars()['ax0'+str(i)].set_ylabel('q', color='blue')
					vars()['ax1'+str(i)].set_ylabel('q', color='blue')
					vars()['ax0'+str(i)].tick_params(axis='y', labelcolor='blue')
					vars()['ax1'+str(i)].tick_params(axis='y', labelcolor='blue')

			else:
	
				labs = [l.get_label() for l in lns]

				vars()['ax'+str(i)].legend(lns, labs, loc='best')
				
				vars()['ax'+str(i)].grid()
				vars()['ax'+str(i)].tick_params(axis='y', labelcolor='red')
				vars()['ax'+str(i)].set_ylabel(r'au', color='red')
				fig.suptitle('Variable %i'%(i+1))
				
				if eq != None:
					vars()['_ax'+str(i)] = vars()['ax'+str(i)].twinx()

					vars()['_ax'+str(i)].plot(self.grid.S,eq.q, '-', color='blue', lw=2, label='q')
					vars()['_ax'+str(i)].set_ylabel('q', color='blue')
					vars()['_ax'+str(i)].tick_params(axis='y', labelcolor='blue')


		#--------------------------------------------------------------------------------
		if show: plt.show()


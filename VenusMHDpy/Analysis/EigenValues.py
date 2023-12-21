import numpy as np
import matplotlib.pyplot as plt
import sys


class EigenValues:
	'''
		Plots and reads eigenvalues from the produced txt file
	'''

	def __init__(self,FileName):
		self.Nradial   = []
		self.Npoloidal = []
		self.ntor      = []
		self.Mmax      = []
		self.Mmin      = []
		self.EV_guess  = []
		self.SimID     = []
		self.EV        = []
		
		
		F = open(FileName,'r')
		lines = F.readlines()
		F.close()

		for i,l in enumerate(lines):
			if 'VENUS-MHDpy' in l:
				j = i+1
				while '#----' not in lines[j]:
					if 'Nradial' in lines[j]:
						self.Nradial.append(int(lines[j].split('=')[1]))
					if 'Npoloidal' in lines[j]:
						self.Npoloidal.append(int(lines[j].split('=')[1]))
					if 'n =' in lines[j]:
						self.ntor.append(int(lines[j].split('=')[1]))
					if 'm =' in lines[j]:
						self.Mmin.append(int(lines[j].split('=')[1].split(',')[0]))
						self.Mmax.append(int(lines[j].split('=')[1].split(',')[1]))
					if 'EV guess' in lines[j]:
						self.EV_guess.append(eval(lines[j].split('=')[1].replace(" ","")))
					if 'Simulation ID' in lines[j]:
						self.SimID.append(lines[j].split(':')[1])
					if 'Eigenvalue list:' in lines[j]:
						temp = []
						k = j+1
						while '#----' not in lines[k]:
							temp.append(eval(lines[k].replace(" ","")))
							k += 1
								
						self.EV.append(temp)
						
					j += 1
		
		self.Nradial   = np.asarray(self.Nradial)
		self.Npoloidal = np.asarray(self.Npoloidal)
		self.ntor      = np.asarray(self.ntor)
		self.Mmax      = np.asarray(self.Mmax)
		self.Mmin      = np.asarray(self.Mmin)
		self.EV_guess  = np.asarray(self.EV_guess)
		self.EV        = np.asarray(self.EV)


if __name__ == '__main__':

	colors = ['k', 'b', 'g', 'y']
	markers = ['x', 'o', '<', '*', '^' ]
	#Open the file
	try:
		for i in range(1,len(sys.argv)):
			j = i-len(colors)*int(i/len(colors))-1
			
			a = EigenValues(sys.argv[i])
			param = np.asarray(a.SimID,dtype=float)
		
			L = len(param)
			if L > 10:
				L = 10
		
			R = np.polyfit(param,a.EV.real,L-1)
			I = np.polyfit(param,a.EV.imag,L-1)
			
			print (np.polyval(R,0.09))
			print (np.polyval(I,0.09))

			plt.figure(1)
			plt.plot(param,a.EV.real,color=colors[j],marker=markers[j], label=sys.argv[i])
			#plt.plot(param,np.polyval(R,param),'--', color='r')
			
			plt.figure(2)
			plt.plot(param,a.EV.imag,color=colors[j],marker=markers[j], label=sys.argv[i])
			#plt.plot(param,np.polyval(I,param),'--', color='r')
			



		plt.figure(1)
		plt.xlabel('Param')
		plt.ylabel(r'$\mathcal{R}$')
		plt.legend(loc='best')
		plt.grid()
		
		plt.figure(2)
		plt.xlabel('Param')
		plt.ylabel(r'$\mathcal{I}$')
		plt.legend(loc='best')
		plt.grid()

		plt.show()
		
	except IndexError:
		None






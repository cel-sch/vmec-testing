import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.insert(0, '../bin/')
import library as lib
from scipy import interpolate


def alphaBspline(x,i,j,X,d=0,sk=10):

	knots = np.concatenate([np.ones(sk,dtype=np.int64)*X[0],X,np.ones(sk,dtype=np.int64)*X[-1]])
	
	R = np.zeros_like(x)
	c = abs(knots[i+sk] - knots[i+sk+j]) >= 1.E-12
	if d == 0: R[c] = (x[c]-knots[i+sk][c])/(knots[i+sk+j][c]-knots[i+sk][c])
	if d == 1: R[c] = 1./(knots[i+sk+j][c]-knots[i+sk][c])

	return R

def Bspline(x,i,j,X,der=0,sk=10):
	f = np.zeros_like(x)

	knots = np.concatenate([np.ones(sk,dtype=np.int64)*X[0],X,np.ones(sk,dtype=np.int64)*X[-1]])

	if j==0:
		if der == 0: f[(knots[i+sk] <= x) & (x < knots[i+sk+1]) | (i == len(X)+j-2) & (x == X[-1])] = 1.
		if der == 1: f = np.zeros_like(x)
	else:
		if der == 0: f = alphaBspline(x,i,j,X)*Bspline(x,i,j-1,X)+(1.-alphaBspline(x,i+1,j,X))*Bspline(x,i+1,j-1,X)
		if der == 1: f = j*(alphaBspline(x,i,j,X,d=1)*Bspline(x,i,j-1,X)-alphaBspline(x,i+1,j,X,d=1)*Bspline(x,i+1,j-1,X))


	return f

def GaussBunching(x,A,r0,sigma):
    fdens = 1.
    
    for i in range(len(r0)):
        
        fdens += A[i]*np.exp(-0.5*((x-r0[i])/sigma[i])**2.)
    
    return fdens

def CreateBunching(Ns,A,r0,sigma, ds_init=0.01):
    
	x_ = np.linspace(0.,1.,1000)
	fdens = GaussBunching(x_,A,r0,sigma)

	fdens0 = Ns/np.trapz(fdens,x_)

	f = [0.]
	ds = ds_init
	while f[-1] < 1.:
		k = 0.5*fdens0*(GaussBunching(f[-1],A,r0,sigma)+GaussBunching(f[-1]+ds,A,r0,sigma))*ds
		ds = ds/k
		
		if f[-1]+ds < 1.0: 
			f.append(f[-1]+ds)
		else:
			if 1.-f[-1] < f[-1]+ds-1.:
				f[-1] = 1.
			else:
				f.append(1.0)

	return np.asarray(f)


if __name__=='__main__':
	'''
		Plots the Bspline basis functions of order j in a non-uniform grid. 
	'''

	j = 4
	der = 1
	N = 6
	#X = np.linspace(0.,1.,N+2)
	X = CreateBunching(N,[5.,3.],[0.2,0.7],[0.05,0.05])
	#X = [0.,0.1,0.2,0.3,0.4,0.42,0.44,0.46,0.48,0.5,0.6,0.7,0.8,0.9,0.91,0.92,0.93,0.94,0.95,1.0]
	x = np.linspace(0.,1.,50000)


	#---------------------------------------------------------------------------
	t0 = time.time()
	
	plt.figure(1)
	for y in X:
	 	plt.axvline(y,color='r')
	
	f  = np.zeros_like(x)
	xvals = []
	for i in range(-j,len(X)-1):
		l = np.ones(len(x),dtype=np.int64)*i
		B1 = Bspline(x,l,j,X)
		f += B1
		plt.plot(x,B1,label=i)
		xvals.append(x[lib.find_nearest(B1,max(abs(B1)))])
		
	#lib.labelLines(plt.gca().get_lines(),fontsize=10,align=False,xvals=xvals)
	
	plt.plot(x,f,'--')
	plt.grid()
	

	
	
	plt.figure(3)
	for y in X:
		plt.axvline(y,color='r')
	
	df = np.zeros_like(x)
	xvals = []
	for i in range(-j,len(X)-1):
		l = np.ones(len(x),dtype=np.int64)*i
		dB1 = Bspline(x,l,j,X,der)
		df += dB1
		plt.plot(x,dB1,label=i)
		xvals.append(x[lib.find_nearest(abs(dB1),max(abs(dB1)))])
	
	# lib.labelLines(plt.gca().get_lines(),fontsize=10,align=False,xvals=xvals)
	
	plt.plot(x,df,'--')
	plt.grid()
	

	
	print ('Time taken: %.4f'%(time.time()-t0))
	#---------------------------------------------------------------------------
	
	
	plt.show()


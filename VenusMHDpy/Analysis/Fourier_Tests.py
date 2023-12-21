import h5py
from scipy.fft import  fft, fftshift, ifft, fftfreq
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, '../../bin')

'''
    Some tests to check the numerical Fourier integration and derivation with scipy fft tool
'''


# Test if the fourier transform works
#==========================================================================================================
C = [-4.,5.,7.,1.,-9.,12.]
S = [6.,-5.,0.2,8.,-4.,11.]
m = [-2.,3.,-4.,5.,0.,7.]


ntheta = 30
T = 2.*np.pi/ntheta
theta  = np.linspace(0.,2.*np.pi,ntheta,endpoint=False) 

A = 0.
dA = 0.
for i in range(len(C)):
    A += (2.+1.j)*C[i]*np.cos(m[i]*theta)+S[i]*np.sin(m[i]*theta)
    dA += -m[i]*(1.0j+2.)*C[i]*np.sin(m[i]*theta)+m[i]*S[i]*np.cos(m[i]*theta)
    
    #IA += (1.0j+2.)*C[i]*np.sin(m[i]*theta)/m[i]-S[i]*np.cos(m[i]*theta)/m[i]

a = fft(A)
f = fftfreq(ntheta,T)

Ap = ifft(2.*np.pi*(1.0j)*f*a)


#plt.plot(theta/(2.*np.pi),dA.real, 'o-')
#plt.plot(theta/(2.*np.pi),dA.imag, '>-')
#plt.plot(theta/(2.*np.pi),Ap.real,'x-')
#plt.plot(theta/(2.*np.pi),Ap.imag,'*-')
#plt.grid()
#plt.show()


#Fourier = fftshift(fftn(f,axes=0), axes=0)/grid.Ntheta
Mmin = -5
Mmax = 5

plt.figure(2)
INT = fftshift(fft(A))/theta.size
m_fft = np.arange( -int(np.fix(theta.size/2)),int(np.ceil(theta.size/2)))


plt.plot(m_fft,INT.imag,'o')
plt.plot(m_fft,INT.real,'x')

        
plt.grid()
plt.show()
#==========================================================================================================


















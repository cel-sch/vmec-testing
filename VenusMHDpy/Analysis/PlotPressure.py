from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import sys
from VenusMHDpy.library import WoutFile,find_nearest, Fsincos, Fsin, Fcos
from scipy import interpolate



wout = WoutFile(sys.argv[1])

#Violates(1) or not(0) stellarator symmetry?
wout.GetVariable('lasym__logical__')

#Poloidal(1) or Toroidal(0) flux as a  radial variable
wout.GetVariable('lrfp__logical__')

#Geometry
#-----------------------------------------------------------
wout.GetVariable('rmnc')
wout.GetVariable('zmns')

if wout.lasym__logical__:
    wout.GetVariable('rmns')
    wout.GetVariable('zmnc')

wout.GetVariable('xm')
wout.GetVariable('xn')
#-----------------------------------------------------------

#Profiles
#-------------------------------------------
wout.GetVariable('ns')
wout.GetVariable('q_factor')
wout.GetVariable('presf')
wout.GetVariable('pres')
wout.GetVariable('am')
wout.GetVariable('at')
wout.GetVariable('ah')
wout.GetVariable('phi')     
wout.GetVariable('chi')
wout.GetVariable('chipf')
wout.GetVariable('jcurv')
wout.GetVariable('jcuru')
wout.GetVariable('bvco')

try:
    wout.GetVariable('machsq')
    M02 = wout.machsq
except KeyError:
    M02 = 0.

s_VMEC = np.linspace(0.,1.,wout.ns)
P           = np.polyval(wout.am[::-1],s_VMEC)
T_norm_VMEC = np.polyval(wout.at[::-1],s_VMEC)
Omega_norm  = np.polyval(wout.ah[::-1],s_VMEC)
U = M02*Omega_norm**2./(2.*wout.rmnc[0][0]**2.*T_norm_VMEC) 


ntheta = 200
theta = np.linspace(0.,2.*np.pi,ntheta)
R = np.zeros((ntheta,wout.ns))
Z = np.zeros((ntheta,wout.ns))

plt.figure()
plt.plot(wout.rmnc[0][0],wout.zmns[0][0],'k+',markersize=10)
for i in range(wout.ns):
    #Reconstruction of the coordinates from Fourier coefficients
    if wout.lasym__logical__:
        R[:,i] = Fsincos(theta,0.,wout.xm,wout.xn,wout.rmnc[i],wout.rmns[i],False)
        Z[:,i] = Fsincos(theta,0.,xm,xn,wout.zmnc[i],wout.zmns[i],False)
    else:
        R[:,i] = Fcos(theta,0.,wout.xm,wout.xn,wout.rmnc[i],False)
        Z[:,i] = Fsin(theta,0.,wout.xm,wout.xn,wout.zmns[i],False)
        
    if i==wout.ns-1 or i%25==0:
        plt.plot(R[:,i],Z[:,i],'k')
        
Pdens = P*np.exp(U*R**2)
Pdens_ = P*np.ones_like(R)

#plt.pcolormesh(R, Z, Pdens, shading='auto')
#plt.contour(R, Z, Pdens_, 30, colors=['red'])
plt.contour(R, Z, Pdens, 25,colors=['blue'], linestyles='dashed')
#plt.contourf(R, Z, Pdens, 3, cmap='hot')
plt.gca().set_aspect('equal')
plt.grid()


plt.figure()
Pbar = wout.presf
plt.plot(s_VMEC,Pbar)
Pbar = np.polyval(wout.am[::-1],s_VMEC)*np.exp(U*wout.rmnc[0][0]**2.)
plt.plot(s_VMEC,Pbar,'--')

plt.grid()


C = np.ones_like(wout.xm)
S = np.zeros_like(wout.xm)
for i,pon in enumerate(wout.xm%2==1):
    if pon:
        C[i] = -1.
        S[i] = 1.
        
P_ = P*np.exp(U*np.sum(wout.rmnc,axis=-1)**2.)
P__ = P[::-1]*np.exp(U[::-1]*np.sum(wout.rmnc*C,axis=-1)[::-1]**2.)

R = np.append(np.sum(wout.rmnc*C,axis=-1)[::-1],np.sum(wout.rmnc,axis=-1))
Z = np.append(-np.sum(wout.zmns*S,axis=-1)[::-1],np.sum(wout.zmns*S,axis=-1))

Pbar = np.append(Pbar[::-1],Pbar)
Prot = np.append(P__,P_)

idx1 = find_nearest(Pbar,max(Pbar))
idx2 = find_nearest(Prot,max(Prot))

plt.figure()
plt.plot(R,Pbar)
plt.axvline(R[idx1])
plt.plot(R,Prot)
plt.axvline(R[idx2])
plt.axvline(wout.rmnc[0][0],linestyle='--',color='r')
plt.grid()

plt.show()

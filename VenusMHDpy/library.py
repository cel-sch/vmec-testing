import numpy as np
from scipy.io import netcdf
from scipy.fft import fftn, ifftn, fftfreq
import matplotlib.pyplot as plt


#*******************************************library.py**********************************************************#
#																												#
'''
	This file contains general routines used in the code.
	
	Author: Guillermo Bustos Ramirez
	Last update: 08/April/2022 

'''

class WoutFile:
	'''
		Class that to open a NetCDF output VMEC file (wout_Label.nc)
	'''
	def __init__(self,filename):
		self.woutfile = filename
	    
	def GetVariable(self,name):
		wout = netcdf.netcdf_file(self.woutfile,'r')
		vars(self)[name] = wout.variables[name][()].copy()
		wout.close()


def find_nearest(array, value,tb=None):
	'''
		find_nearest                                                           
		Function to find the index of the nearest value of an array. It's possible to specify if one wants the 'top'
		or the 'bottom' value (values above or below the specified value)
		                                                                       
		input: array, value, tb                                                    
		output: index of the nearest value in the array to the specified value 
	'''
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	
	if tb==None: return idx
	if tb=='top':
		if array[idx] >= value:
			return idx
		else:
			return idx+1
	if tb=='bottom':
		if array[idx] <= value:
			return idx
		else:
			return idx-1

def PlotParameters(text=999, axestitle=24, xyaxes=22, xtick=16, ytick=16, legend=20, title=20):
    plt.rc('font', size=text)          # controls default text sizes
    plt.rc('axes', titlesize=axestitle)     # fontsize of the axes title
    plt.rc('axes', labelsize=xyaxes)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=xtick)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=ytick)    # fontsize of the tick labels
    plt.rc('legend', fontsize=legend)    # legend fontsize
    #plt.rc('figure', titlesize=title)  # fontsize of the figure title

    return None



#**********************************************************************************************************
#** Fcos, Fsin, dFcos_u, dFsin_u                                                                         **
#** Given the coefficients of the Fourier harmonics of a function or its derivatives (cos coefficients   **
#** for Fcos and sin coefficients for Fsin) the function returns an array of the values of the function  **
#** of which harmonics were provided.                                                                    ** 
#**                                                                                                      **
#** input  u,v  -> arrays of theta and phi (u,v respectivelly) of the fourier harmonics                  **
#**        m,n  -> arrays of the poloidal (m) and toroidal (n) modes (Given by VMEC)                     **
#**        coeff  -> array with the fourier coefficients of the desired function (Given by VMEC)         **
#**                  Enter rmns[s] for example, where 's' is the index of the flux surface               **
#**        flag  -> Flag to set 3D (True) or 2D (False) equilibrium. 3D uses all mode numbers,           **
#**                 while 2D uses only n=0 harmonics.                                                    **
#** output: array with the function as a function of theta and phi                                       **
#**********************************************************************************************************
def Fcos(u,v,m,n,coeff,flag):
    if flag == False:
        F = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            F += coeff[i]*np.cos(m[i]*u)
        return F
    else:
        try:
            F = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    F[j] += coeff[i]*np.cos(m[i]*u-n[i]*v[j])
        except TypeError:
            F = 0
            for i in range(len(m)):
                F += coeff[i]*np.cos(m[i]*u-n[i]*v)
        return np.asarray(F)
    
def Fsin(u,v,m,n,coeff,flag):
    if flag == False:
        F = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            F += coeff[i]*np.sin(m[i]*u)
        return F
    else:
        try:
            F = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    F[j] += coeff[i]*np.sin(m[i]*u-n[i]*v[j])
        except TypeError:
            F = 0
            for i in range(len(m)):
                F += coeff[i]*np.sin(m[i]*u-n[i]*v)
        return np.asarray(F)

def dFcos_u(u,v,m,n,coeff,flag):
    if flag == False:
        dF = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            dF += coeff[i]*m[i]*np.sin(m[i]*u)
        return -dF
    else:
        try:
            dF = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    dF[j] += coeff[i]*m[i]*np.sin(m[i]*u-n[i]*v[j])
        except TypeError:
            dF = 0
            for i in range(len(m)):
                dF += coeff[i]*m[i]*np.sin(m[i]*u-n[i]*v)
        return -np.asarray(dF)
    
def dFcos_v(u,v,m,n,coeff,flag):
    if flag == False:
        dF = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            dF += 0.0
        return -dF
    else:
        try:
            dF = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    dF[j] += coeff[i]*n[i]*np.sin(m[i]*u-n[i]*v[j])
        except TypeError:
            dF = 0
            for i in range(len(m)):
                dF += coeff[i]*n[i]*np.sin(m[i]*u-n[i]*v)
        return np.asarray(dF)

def dFsin_u(u,v,m,n,coeff,flag):
    if flag == False:
        dF = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            dF += coeff[i]*m[i]*np.cos(m[i]*u)
        return dF
    else:
        try:
            dF = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    dF[j] += coeff[i]*m[i]*np.cos(m[i]*u-n[i]*v[j])
        except TypeError:
            dF = 0
            for i in range(len(m)):
                dF += coeff[i]*m[i]*np.cos(m[i]*u-n[i]*v)
        return np.asarray(dF)
    
def dFsin_v(u,v,m,n,coeff,flag):
    if flag == False:
        dF = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            dF += 0.
        return dF
    else:
        try:
            dF = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    dF[j] += coeff[i]*n[i]*np.cos(m[i]*u-n[i]*v[j])
        except TypeError:
            dF = 0
            for i in range(len(m)):
                dF += coeff[i]*n[i]*np.cos(m[i]*u-n[i]*v)
        return -np.asarray(dF)
    
def Fsincos(u,v,m,n,C,S,flag):
    if flag == False:
        F = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            F += C[i]*np.cos(m[i]*u)+S[i]*np.sin(m[i]*u)
        return F
    else:
        try:
            F = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    F[j] += C[i]*np.cos(m[i]*u-n[i]*v[j])+S[i]*np.sin(m[i]*u-n[i]*v[j])
        except TypeError:
            F = 0
            for i in range(len(m)):
                F += C[i]*np.cos(m[i]*u-n[i]*v)+S[i]*np.sin(m[i]*u-n[i]*v)
        return np.asarray(F)

def dFsincos_u(u,v,m,n,C,S,flag):
    if flag == False:
        F = 0
        index0 = np.where(n == 0)[0]
        for i in index0:
            F += -C[i]*m[i]*np.sin(m[i]*u)+S[i]*m[i]*np.cos(m[i]*u)
        return F
    else:
        try:
            F = [0]*len(v)
            for j in range(len(v)):
                for i in range(len(m)):
                    F[j] += -C[i]*m[i]*np.sin(m[i]*u-n[i]*v[j])+S[i]*m[i]*np.cos(m[i]*u-n[i]*v[j])
        except TypeError:
            F = 0
            for i in range(len(m)):
                F += -C[i]*m[i]*np.sin(m[i]*u-n[i]*v)+S[i]*m[i]*np.cos(m[i]*u-n[i]*v)
        return np.asarray(F)


#****************************************************************************************************************
#**Fourier                                                                                                     **
#**                                                                                                            **
#** Fourier decompose a 1D funcion and returns an array with the amplitudes of the Fourier coefficients        **
#** in the order (-fn,fn). So to reconstruct the function, one has to loop as:                                 **
#** for i in range(int(t.size/2)-1-Mmax,int(theta.size/2)-1+Mmax):                                             **
#**     Reconstruction += np.real(coeff[i])*np.cos(frequency[i]*t)-np.imag(coeff[i])*np.sin(frequency[i]*t)    **
#**                                                                                                            **
#** Input:  array --> Array of values to Fourier analyze                                                       **
#** Output:  A/t.size --> Amplitude of the fourier  coefficients                                               **
#**          frequencies --> Corresponding frequencies                                                         **
#****************************************************************************************************************
def Fourier(array, t):
    
    f = np.fft.fft(array)
    if t.size%2==0:
        A = np.zeros(t.size-1,dtype=complex)
        A[0:int(t.size/2)-1] = f[int(t.size/2)+1:]
        A[int(t.size/2)-1] = f[0]
        A[int(t.size/2):] = f[1:int(t.size/2)]
        freq = np.arange(-int(t.size/2)+1,int(t.size/2),1)

    else:
        A = np.zeros(t.size,dtype=complex)
        A[0:(t.size-1)//2] = f[(t.size-1)//2+1:]
        A[t.size//2] = f[0]
        A[(t.size+1)//2:] = f[1:(t.size+1)//2]
        freq = np.arange(-(t.size-1)//2,(t.size+1)//2,1)
   
    return A/t.size, freq


#****************************************************************************************************************
#**GaussQuadrature                                                                                             **
#**                                                                                                            **
#** Stores the coefficients to calculate the position of the Gauss nodes, as well as the weigths for each node **
#**                                                                                                            **
#****************************************************************************************************************

class GaussQuadrature:
    def __init__(self,p):
        
        #Select the positiosn of the Gauss nodes.
        if p == 1: self.Pcoeff = [0.]
        if p == 2: self.Pcoeff = [-np.sqrt(1./3.),np.sqrt(1./3.)]
        if p == 3: self.Pcoeff = [-np.sqrt(3./5.),0.,np.sqrt(3./5.)]
        if p == 4: self.Pcoeff = [-np.sqrt(3./7.+2.*np.sqrt(6./5.)/7.),-np.sqrt(3./7.-2.*np.sqrt(6./5.)/7.),np.sqrt(3./7.-2.*np.sqrt(6./5.)/7.),np.sqrt(3./7.+2.*np.sqrt(6./5.)/7.)]
        if p == 5: self.Pcoeff = [-np.sqrt(5.+2.*np.sqrt(10./7))/3., -np.sqrt(5.-2.*np.sqrt(10./7))/3., 0., np.sqrt(5.-2.*np.sqrt(10./7))/3., np.sqrt(5.+2.*np.sqrt(10./7))/3.]
    
        self.Weight = np.ones(shape=(p,1))	
        if p == 1: self.Weight[:,0] = 2.
        if p == 2: self.Weight[:,0] = [1.,1.]
        if p == 3: self.Weight[:,0] = [5./9.,8./9.,5./9.]
        if p == 4: self.Weight[:,0] = [(18.-np.sqrt(30.))/36.,(18.+np.sqrt(30.))/36.,(18.+np.sqrt(30.))/36.,(18.-np.sqrt(30.))/36.]
        if p == 5: self.Weight[:,0] = [(322.-13.*np.sqrt(70.))/900.,(322.+13.*np.sqrt(70.))/900., 128./225., (322.+13.*np.sqrt(70.))/900.,(322.-13.*np.sqrt(70.))/900.]
       

 
class GaussJacobi:
    def __init__(self,p):
        
        #Select the positiosn of the Gauss nodes.
        if p == 2: self.Pcoeff = [-0.999995000006250,0.333337222198843]
        if p == 3: self.Pcoeff = [-0.999997777775309,-0.289894139677351,0.689899117438264]
        if p == 4: self.Pcoeff = [-0.999998749997266,-0.575316121535975,0.181068260268316,0.822824570439263]
        if p == 5: self.Pcoeff = [-0.999999199997760,-0.720478227479051,-0.167178894832870,0.446315083462093,0.885791856124685]
    
        self.Weight = np.ones(shape=(p,1))	
        if p == 2: self.Weight[:,0] = [99999.5682,1.12499420424522]
        if p == 3: self.Weight[:,0] = [99998.804,1.44339981684823,0.44547474733231]
        if p == 4: self.Weight[:,0] = [99998.245,1.54864169207012,0.65735975042087,0.24189156392928]
        if p == 5: self.Weight[:,0] = [99997.807,1.5963072438726,0.74884187938281,0.38906675121166,0.15241784520583]
        
#****************************************************************************************************************
#**SplineDer                                                                                                   **
#**                                                                                                            **
#** Takes the derivative of a function with respect to s or theta (or both) using spline representation with   **
#** respect to s, and a fourier derivative with respect to theta.                                              **
#**                                                                                                            **
#****************************************************************************************************************

def SplineDer1(fun,grid,der_s):

	result = fun.copy()
	for d in range(der_s):
		result = np.gradient(result,grid.S,edge_order=2)
	
	return result

def SplineDer2(fun,grid,der_s,der_theta):

	result = fun.copy()
	for d in range(der_s):
		result = np.gradient(result,grid.S,axis=1,edge_order=2)
		
	if der_theta != 0:
		T = fftn(result,axes=0)
		freq = fftfreq(grid.Ntheta,2.*np.pi/grid.Ntheta)
		freq = np.reshape(freq,(freq.size,1))
		result = ifftn(2.*np.pi*(1.0j)*freq*T,axes=0)
        
	return result.real


#****************************************************************************************************************
#**labelLine                                                                                                   **
#**                                                                                                            **
#** Label lines more easily direclty on the line, and not in a square. This was obtainedd from stack overflow **
#** https://stackoverflow.com/questions/16992038/inline-labels-in-matplotlib                                   **
#**                                                                                                            **
#****************************************************************************************************************


#Label line with line2D label data
def labelLine(line,x,label=None,align=True,**kwargs):

    ax = line.axes
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy,dx))

        #Transform to screen co-ordinates
        pt = np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_facecolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x,y,label,rotation=trans_angle,**kwargs)

def labelLines(lines,align=True,xvals=None,**kwargs):

    ax = lines[0].axes
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin,xmax = ax.get_xlim()
        xvals = np.linspace(xmin,xmax,len(labLines)+2)[1:-1]

    for line,x,label in zip(labLines,xvals,labels):
        labelLine(line,x,label,align,**kwargs)



def has_twin(ax):
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            return True
    return False


        
        



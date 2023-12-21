from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../bin/')
import time
import library as lib
from scipy import interpolate



def Plot(Files):
    '''
        Plots the main equilibrium quantities of a VMEC wout_*.nc file. It can be used to compare
        different files, and can be called from terminal.
    '''
    
    
    
    #Read the output file
    #========================================================================================
    colors = ['k', 'b', 'r', 'g', 'y']

    print ('Reading files...')

    for i in range(len(Files)):
        vars()['F'+str(i)]= lib.WoutFile(Files[i])
        
        #Violates(1) or not(0) stellarator symmetry?
        vars()['F'+str(i)].GetVariable('lasym__logical__')
    
        #Poloidal(1) or Toroidal(0) flux as a  radial variable
        vars()['F'+str(i)].GetVariable('lrfp__logical__')

        #Geometry
        #-----------------------------------------------------------
        vars()['F'+str(i)].GetVariable('rmnc')
        vars()['F'+str(i)].GetVariable('zmns')
        
        if vars()['F'+str(i)].lasym__logical__:
            vars()['F'+str(i)].GetVariable('rmns')
            vars()['F'+str(i)].GetVariable('zmnc')
        
        vars()['F'+str(i)].GetVariable('xm')
        vars()['F'+str(i)].GetVariable('xn')
        #-----------------------------------------------------------
        
        #Profiles
        #-------------------------------------------
        vars()['F'+str(i)].GetVariable('ns')
        vars()['F'+str(i)].GetVariable('q_factor')
        vars()['F'+str(i)].GetVariable('presf')
        vars()['F'+str(i)].GetVariable('pres')
        vars()['F'+str(i)].GetVariable('am')
        vars()['F'+str(i)].GetVariable('at')
        vars()['F'+str(i)].GetVariable('ah')
        vars()['F'+str(i)].GetVariable('phi')     
        vars()['F'+str(i)].GetVariable('chi')
        vars()['F'+str(i)].GetVariable('chipf')
        vars()['F'+str(i)].GetVariable('jcurv')
        vars()['F'+str(i)].GetVariable('jcuru')
        vars()['F'+str(i)].GetVariable('bvco')
        
        try:
            vars()['F'+str(i)].GetVariable('machsq')
            M02 = vars()['F'+str(i)].machsq
        except KeyError:
            M02 = 0.
        
        #Set radial variable for plots
        setRadial = 'SqrtNormPolFlux'
        if setRadial == 'NormTorFlux':
            S = vars()['F'+str(i)].phi/vars()['F'+str(i)].phi[-1]
            rlabel = r'$\Psi_T/\Psi_T$'
            f = -vars()['F'+str(i)].phi[-1]/vars()['F'+str(i)].q_factor
        if setRadial == 'SqrtNormTorFlux':
            S = np.sqrt(vars()['F'+str(i)].phi/vars()['F'+str(i)].phi[-1])
            rlabel = r'$\sqrt{\Psi_T/\Psi_T}$'
            f = -2.*S*vars()['F'+str(i)].phi[-1]/vars()['F'+str(i)].q_factor
        if setRadial == 'NormPolFlux':
            S = vars()['F'+str(i)].chi/vars()['F'+str(i)].chi[-1]
            rlabel = r'$\Psi_P/\Psi_P$'
            f = -vars()['F'+str(i)].chi[-1]*np.ones_like(S)
        if setRadial == 'SqrtNormPolFlux':
            S = np.sqrt(vars()['F'+str(i)].chi/vars()['F'+str(i)].chi[-1])
            rlabel = r'$\sqrt{\Psi_P/\Psi_P}$'
            f = -2.*S*vars()['F'+str(i)].chi[-1]
        #-------------------------------------------

        #Other parameters
        #--------------------------------------
        vars()['F'+str(i)].GetVariable('b0')
        vars()['F'+str(i)].GetVariable('Rmajor_p')
        #--------------------------------------
        

        #Set parameters on the plots.
        lib.PlotParameters(text=12, axestitle=12, xyaxes=12, xtick=10, ytick=10, legend=12, title=12)

        j = i-len(colors)*int(i/len(colors))    

        #Boundary fixed poloidal plane
        #--------------------------------------------------------------------
        Xsize = 200
        X = np.linspace(0.,2.*np.pi, Xsize)

        #2D surfaces
        if vars()['F'+str(i)].lasym__logical__:
            R = lib.Fsincos(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].rmnc[-1],vars()['F'+str(i)].rmns[-1],False)
            Z = lib.Fsincos(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].zmnc[-1],vars()['F'+str(i)].zmns[-1],False)
            
            R0 = lib.Fsincos(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].rmnc[0],vars()['F'+str(i)].rmns[0],False)
            Z0 = lib.Fsincos(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].zmnc[0],vars()['F'+str(i)].zmns[0],False)
        else:	
            R = lib.Fcos(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].rmnc[-1],False)
            Z = lib.Fsin(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].zmns[-1],False)

            R0 = lib.Fcos(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].rmnc[0],False)
            Z0 = lib.Fsin(X,0.,vars()['F'+str(i)].xm,vars()['F'+str(i)].xn,vars()['F'+str(i)].zmns[0],False)
            #Plot the LCFS 2D at the phi angle
            
        plt.figure(1)
        plt.plot(R,Z, color=colors[j], label=Files[i])
        plt.plot(R0,Z0,'+', color=colors[j])
        


        plt.title(r'2D LCFS')
        plt.xlabel(r'R [m]')
        plt.ylabel(r'Z [m]')
        plt.gca().set_aspect('equal')
        #--------------------------------------------------------------------
        
        #Calulate the pressure, temperature and rotation profiles from polynomial coefficients
        #------------------------------------------------------
        s_VMEC = np.linspace(0.,1.,vars()['F'+str(i)].ns)
        P          = np.polyval(vars()['F'+str(i)].am[::-1],s_VMEC)
        T_norm     = np.polyval(vars()['F'+str(i)].at[::-1],s_VMEC)
        Omega_norm = np.polyval(vars()['F'+str(i)].ah[::-1],s_VMEC)
        
        U = M02*Omega_norm**2./(2.*vars()['F'+str(i)].rmnc[0][0]**2.*T_norm) 
        #------------------------------------------------------
        
        
        #Calculate the ballooning alpha parameter and shear
        #------------------------------------------------------
        NormTorFlux = vars()['F'+str(i)].phi/vars()['F'+str(i)].phi[-1]
        pres_spline = interpolate.splrep(np.sqrt(NormTorFlux),P)
        dpresf = interpolate.splev(np.sqrt(NormTorFlux),pres_spline,der=1) 
        mu0 = 4.*np.pi*10.**(-7.)
        alpha = -2.*mu0*dpresf*vars()['F'+str(i)].q_factor**2.*vars()['F'+str(i)].rmnc[0][0]/vars()['F'+str(i)].b0**2.
        beta  = 2.*mu0*P/vars()['F'+str(i)].b0**2.


        print ("%s B0 = %.4f"%(Files[i],vars()['F'+str(i)].b0))
        print ("%s R0 = %.4f"%(Files[i],vars()['F'+str(i)].rmnc[0][0]))


        q_spline = interpolate.splrep(S,vars()['F'+str(i)].q_factor)
        dq = interpolate.splev(S,q_spline,der=1)

        shear = S*dq/vars()['F'+str(i)].q_factor
        #------------------------------------------------------
        
        

        plt.figure(2)
        plt.plot(S,vars()['F'+str(i)].q_factor, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel('q')
        
        plt.figure(3)
        plt.plot(S,P, label=Files[i], color=colors[j])
        plt.plot(S,vars()['F'+str(i)].presf,'--', label=Files[i], color='y')
        #plt.plot(S,vars()['F'+str(i)].presf*np.exp(-U*vars()['F'+str(i)].rmnc[0][0]**2.),'--', label=Files[i], color='y')
        plt.xlabel(rlabel)
        plt.ylabel('P')
        
        plt.figure(4)
        plt.plot(S,vars()['F'+str(i)].jcuru, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$<j_u>$')
        
        plt.figure(5)
        plt.plot(S,vars()['F'+str(i)].jcurv, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$<j_v>$')
        
        plt.figure(6)
        plt.plot(S,alpha, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$\alpha$')
        
        plt.figure(7)
        plt.plot(S,shear, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$s$')
        
        plt.figure(8)
        plt.plot(S[1:],vars()['F'+str(i)].bvco[1:], label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$F$')
        
        plt.figure(9)
        plt.plot(S,vars()['F'+str(i)].phi, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$\Psi_T$')
        
        plt.figure(10)
        plt.plot(S,vars()['F'+str(i)].chi, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$\Psi_P$')
        
        plt.figure(11)
        plt.plot(S,f/(2.*np.pi), label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$f$')
        
        plt.figure(12)
        plt.plot(S,T_norm, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$\hat{T}$')
        
        plt.figure(13)
        plt.plot(S,Omega_norm, label=Files[i], color=colors[j])
        plt.xlabel(rlabel)
        plt.ylabel(r'$\hat{\Omega}$')

    for i in range(1,14):
        plt.figure(i)
        plt.grid()
        plt.ticklabel_format(axis="y", style="sci",scilimits=(0,0))
        plt.legend()


    plt.show()


if __name__ == '__main__':

	#Open the file
	try:
		woutfiles = []
		for i in range(1,len(sys.argv)):
			woutfiles.append(sys.argv[i])

		Plot(woutfiles)

	except IndexError:
		None






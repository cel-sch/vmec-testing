import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.alpha*eq.dm) - eq.alpha*np.cos(eq.r/(1 + eq.delta - eq.r**2)) + eq.alpha*eq.shear*np.sin(eq.r/(1 + eq.delta - eq.r**2)) + (eq.alpha*eq.r*eq.shear*np.sin(eq.r/(1 + eq.delta - eq.r**2)))/(-1 - eq.delta + eq.r**2) + eq.alpha**2*np.sin(eq.r/(1 + eq.delta - eq.r**2))**2



        


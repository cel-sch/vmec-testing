import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.alpha*eq.dm) - eq.alpha*np.cos(eq.r) + eq.alpha*eq.shear*np.sin(eq.r) - eq.alpha*eq.r*eq.shear*np.sin(eq.r) + eq.alpha**2*np.sin(eq.r)**2



        


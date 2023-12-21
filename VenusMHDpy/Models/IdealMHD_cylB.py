import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = eq.rho0/2.



        

class M_01:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_02:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_10:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_11:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.rho0*eq.S**2)/2.



        

class M_12:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz*eq.rho0*eq.S**2)/(2.*eq.q)



        

class M_20:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_21:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz*eq.rho0*eq.S**2)/(2.*eq.q)



        

class M_22:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz**2*eq.rho0)/2. + (eq.Bz**2*eq.rho0*eq.S**2)/(2.*eq.q**2)



        


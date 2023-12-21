import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.F*eq.g12**2*eq.rho)/(eq.B2*eq.g*eq.q**3*eq.R2)) + (eq.g11*eq.R2*eq.rho)/(eq.F*eq.g*eq.q)



        

class M_01:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.F*eq.g12*eq.rho*imag)/(eq.B2*eq.g*eq.q))



        

class M_02:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.g12*eq.rho*imag)/(eq.g*eq.q))



        

class M_10:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.F*eq.g12*eq.rho*imag)/(eq.B2*eq.g*eq.q)



        

class M_11:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.F*eq.q*eq.R2*eq.rho)/(eq.B2*eq.g)) + (eq.q*eq.R2**2*eq.rho)/(eq.F*eq.g)



        

class M_12:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.q*eq.R2*eq.rho)/eq.g) + (eq.B2*eq.q*eq.R2**2*eq.rho)/(eq.F**2*eq.g)



        

class M_20:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_21:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_22:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.B2*eq.q*eq.R2*eq.rho)/(eq.F*eq.g)



        


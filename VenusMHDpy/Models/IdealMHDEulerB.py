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



        

class M_03:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_04:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_05:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_06:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_07:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

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



        

class M_13:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_14:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_15:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_16:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_17:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

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



        

class M_23:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_24:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_25:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_26:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_27:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_30:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_31:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_32:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_33:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = 1



        

class M_34:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_35:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_36:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_37:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_40:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_41:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_42:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_43:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_44:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = 1



        

class M_45:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_46:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_47:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_50:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_51:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_52:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_53:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_54:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_55:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = 1



        

class M_56:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_57:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_60:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_61:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_62:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_63:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_64:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_65:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_66:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.q*eq.R2)/(eq.F*eq.g)



        

class M_67:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_70:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_71:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_72:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_73:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_74:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_75:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_76:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_77:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.q*eq.R2*eq.rho)/(eq.F*eq.g*(-1 + eq.gamma))



        


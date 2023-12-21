import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_01:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_02:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_03:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.g11*eq.R2*eq.rhorot)/(eq.F*eq.g*eq.q)



        

class M_04:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.g12*eq.R2*eq.rhorot*imag)/(eq.F*eq.g*eq.q))



        

class M_05:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.g12*eq.rhorot)/(eq.g*eq.q)



        

class M_10:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_11:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_12:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_13:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.g12*eq.R2*eq.rhorot*imag)/(eq.F*eq.g*eq.q)



        

class M_14:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.q*eq.R2**2*eq.rhorot)/(eq.F*eq.g)) + (eq.B2*eq.q*eq.R2**3*eq.rhorot)/(eq.F**3*eq.g)



        

class M_15:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (-((eq.q*eq.R2*eq.rhorot)/eq.g) + (eq.B2*eq.q*eq.R2**2*eq.rhorot)/(eq.F**2*eq.g))*imag



        

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
        self.DerBs = [False,False,False,False]




        

class M_23:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.g12*eq.rhorot)/(eq.g*eq.q)



        

class M_24:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = ((eq.q*eq.R2*eq.rhorot)/eq.g - (eq.B2*eq.q*eq.R2**2*eq.rhorot)/(eq.F**2*eq.g))*imag



        

class M_25:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.B2*eq.q*eq.R2*eq.rhorot)/(eq.F*eq.g)



        

class M_30:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,True,True]

        self.DerPol_DerBs00 = [[0,0,0,1,1,],[0,0,1,0,1,],[-3,-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(5,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.c2*eq.dg11du*eq.dR2du*eq.F*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q**3) + (eq.c2*eq.dg12ds*eq.dR2du*eq.F*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**3) + (eq.c2*eq.dg11du*eq.dR2ds*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q**3) - (eq.c2*eq.dg12ds*eq.dR2ds*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**3) - (eq.c2*eq.dg11ds*eq.dR2du*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q**3) + (eq.c2*eq.dg11ds*eq.dR2ds*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q**3)
        self.var_D00[1] = (eq.dqds**2*eq.F)/(eq.g*eq.q**3) + (eq.dFds*eq.dqds)/(eq.g*eq.q**2) - (eq.dProtds*eq.dR2ds)/(eq.F*eq.g*eq.q) + (2*eq.dFds*eq.dR2ds*eq.gamma*eq.Prot)/(eq.F**2*eq.g*eq.q) - (eq.dR2ds**2*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q*eq.R2) - (eq.B2*eq.dqds**2*eq.R2)/(eq.F*eq.g*eq.q**3) + (eq.dProtds*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) + (eq.dFds*eq.dProtds*eq.R2)/(eq.F**2*eq.g*eq.q) - (eq.c2*eq.dR2ds*eq.drhorotds*eq.Omega**2*eq.R2)/(2.*eq.F*eq.g*eq.q) - (eq.dFds**2*eq.gamma*eq.Prot*eq.R2)/(eq.F**3*eq.g*eq.q) + (eq.c0*eq.dR2ds**2*eq.Omega**2*eq.rhorot)/(4.*eq.F*eq.g*eq.q) - (eq.c1*eq.dR2ds**2*eq.Omega**2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) - (eq.c2*eq.dOmegads*eq.dR2ds*eq.Omega*eq.R2*eq.rhorot)/(eq.F*eq.g*eq.q) + (eq.c1*eq.dFds*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F**2*eq.g*eq.q) - (eq.c2*eq.d2R2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) + (eq.c2*eq.dg11ds*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.F*eq.g*eq.g11*eq.q) - (eq.dqds*eq.dR2ds*eq.Prot*eq.R2*eq.U)/(eq.F*eq.g*eq.q**2) - (eq.F*eq.g11*ntor**2)/(eq.g*eq.q*eq.R2) + (eq.c0*eq.g11*eq.Omega**2*eq.R2*eq.rhorot*ntor**2)/(eq.F*eq.g*eq.q)
        self.var_D00[2] = -((eq.dqds*eq.F*eq.g12)/(eq.g*eq.q**4*eq.R2)) + (eq.F*eq.g11*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[3] = -((eq.dqds*eq.F*eq.g12)/(eq.g*eq.q**4*eq.R2)) - (eq.F*eq.g11*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[4] = -((eq.F*eq.g11)/(eq.g*eq.q**3*eq.R2))

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.dqds*eq.F)/(eq.g*eq.q**2)) - (eq.dR2ds*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q) + (eq.B2*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) - (eq.dProtds*eq.R2)/(eq.F*eq.g*eq.q) + (eq.dFds*eq.gamma*eq.Prot*eq.R2)/(eq.F**2*eq.g*eq.q) - (eq.c1*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) + (eq.dR2ds*eq.Prot*eq.R2*eq.U)/(eq.F*eq.g*eq.q) - (eq.F*eq.g12*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D01[1] = (eq.F*eq.g12)/(eq.g*eq.q**3*eq.R2)

        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.dqds*eq.F)/(eq.g*eq.q**2)) - (eq.dR2ds*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q) + (eq.B2*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) - (eq.dProtds*eq.R2)/(eq.F*eq.g*eq.q) + (eq.dFds*eq.gamma*eq.Prot*eq.R2)/(eq.F**2*eq.g*eq.q) + (eq.F*eq.g12*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D10[1] = (eq.F*eq.g12)/(eq.g*eq.q**3*eq.R2)

        self.DerPol_DerBs11 = [[0,],[0,],[-1,]]

        self.var_D11 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D11[0] = -((eq.B2*eq.R2)/(eq.F*eq.g*eq.q)) - (eq.gamma*eq.Prot*eq.R2)/(eq.F*eq.g*eq.q)
        

class M_31:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,0,0,1,],[0,0,1,0,],[-3,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = ((eq.c2*eq.dg11du*eq.dR2du*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q**3) - (eq.c2*eq.dg11du*eq.dR2ds*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q**3) + (eq.c2*eq.dR2ds*eq.dR2du*eq.F*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.c2*eq.dR2ds**2*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.c2*eq.dqds*eq.dR2du*eq.F*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) - (eq.c2*eq.dqds*eq.dR2ds*eq.F*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dR2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dR2ds**2*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q**2) + (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q**2) + (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g**3*eq.q) - (eq.c2*eq.dB2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g**3*eq.q) + (eq.c2*eq.dB2ds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.q))*imag
        self.var_D00[1] = (eq.dFds*ntor)/eq.g + (eq.dqds*eq.F*ntor)/(eq.g*eq.q) + (eq.dProtds*eq.R2*ntor)/(eq.F*eq.g) - (eq.B2*eq.dqds*eq.R2*ntor)/(eq.F*eq.g*eq.q) - (eq.dR2ds*eq.Prot*eq.R2*eq.U*ntor)/(eq.F*eq.g) + imag*((eq.dProtdu*eq.dR2ds)/(eq.F*eq.g*eq.q) - (eq.dFds*eq.dR2du*eq.gamma*eq.Prot)/(eq.F**2*eq.g*eq.q) + (eq.dR2ds*eq.dR2du*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q*eq.R2) - (eq.dFds*eq.dProtdu*eq.R2)/(eq.F**2*eq.g*eq.q) + (eq.c2*eq.dR2ds*eq.drhorotdu*eq.Omega**2*eq.R2)/(2.*eq.F*eq.g*eq.q) - (eq.c0*eq.dR2ds*eq.dR2du*eq.Omega**2*eq.rhorot)/(4.*eq.F*eq.g*eq.q) + (eq.c1*eq.dR2ds*eq.dR2du*eq.Omega**2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) + (eq.c2*eq.d2R2dsdu*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) - (eq.c2*eq.dg11du*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.F*eq.g*eq.g11*eq.q) + (eq.F*eq.g12*ntor**2)/(eq.g*eq.q*eq.R2) - (eq.c0*eq.g12*eq.Omega**2*eq.R2*eq.rhorot*ntor**2)/(eq.F*eq.g*eq.q))
        self.var_D00[2] = (-(eq.dFds/(eq.g*eq.q)) + (eq.dR2ds*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q) - (eq.dFds*eq.gamma*eq.Prot*eq.R2)/(eq.F**2*eq.g*eq.q) + (eq.c1*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q))*imag
        self.var_D00[3] = -((eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2))


        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = ((eq.dR2du*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q) + (eq.dProtdu*eq.R2)/(eq.F*eq.g*eq.q))*imag - (eq.F*ntor)/eq.g + (eq.B2*eq.R2*ntor)/(eq.F*eq.g)
        self.var_D10[1] = (eq.F/(eq.g*eq.q) + (eq.gamma*eq.Prot*eq.R2)/(eq.F*eq.g*eq.q))*imag

        

class M_32:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,0,0,],[0,0,1,],[-3,-1,-1,]]

        self.var_D00 = np.zeros(shape=(3,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.c2*eq.dqds*eq.dR2du*eq.F**2*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.c2*eq.dqds*eq.dR2ds*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.B2*eq.c2*eq.dR2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dR2ds**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q) - (eq.c2*eq.dg11du*eq.dR2du*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q**3*eq.R2) + (eq.c2*eq.dg11du*eq.dR2ds*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q**3*eq.R2) - (eq.c2*eq.dR2ds*eq.dR2du*eq.F**2*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q*eq.R2) + (eq.c2*eq.dR2ds**2*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q*eq.R2) + (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.c2*eq.dB2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q)
        self.var_D00[1] = (eq.dFds*eq.dProtdu)/(eq.F*eq.g*eq.q) - (eq.c2*eq.dR2ds*eq.drhorotdu*eq.Omega**2)/(2.*eq.g*eq.q) - (eq.dProtdu*eq.dR2ds)/(eq.g*eq.q*eq.R2) - (eq.c2*eq.d2R2dsdu*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q) + (eq.c2*eq.dg11du*eq.dR2ds*eq.Omega**2*eq.rhorot)/(4.*eq.g*eq.g11*eq.q) + (eq.c0*eq.dR2ds*eq.dR2du*eq.Omega**2*eq.rhorot)/(4.*eq.g*eq.q*eq.R2) + (eq.c0*eq.g12*eq.Omega**2*eq.rhorot*ntor**2)/(eq.g*eq.q) + imag*((eq.dFds*eq.gamma*eq.Prot*ntor)/(eq.F*eq.g) - (eq.dR2ds*eq.gamma*eq.Prot*ntor)/(eq.g*eq.R2) + (eq.c0*eq.dR2ds*eq.Omega**2*eq.rhorot*ntor)/eq.g - (eq.c1*eq.dR2ds*eq.Omega**2*eq.rhorot*ntor)/(2.*eq.g))
        self.var_D00[2] = (eq.dFds*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q) - (eq.dR2ds*eq.gamma*eq.Prot)/(eq.g*eq.q*eq.R2) - (eq.c1*eq.dR2ds*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q)


        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -(eq.dProtdu/(eq.g*eq.q)) - (eq.gamma*eq.Prot*imag*ntor)/eq.g
        self.var_D10[1] = -((eq.gamma*eq.Prot)/(eq.g*eq.q))

        

class M_33:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.F*eq.g12**2*eq.kappa*eq.Kpar*eq.rhorot*eq.Uthi)/(eq.B2*eq.g*eq.q**3*eq.R*eq.R2)) - (2*eq.c3*eq.g11*eq.Omega*eq.R2*eq.rhorot*imag*ntor)/(eq.F*eq.g*eq.q)



        

class M_34:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (-((eq.F*eq.g12*eq.kappa*eq.Kpar*eq.rhorot*eq.Uthi)/(eq.B2*eq.g*eq.q*eq.R)) + (eq.g12*eq.kappa*eq.Kpar*eq.R2*eq.rhorot*eq.Uthi)/(eq.F*eq.g*eq.q*eq.R))*imag - (2*eq.c3*eq.g12*eq.Omega*eq.R2*eq.rhorot*ntor)/(eq.F*eq.g*eq.q)



        

class M_35:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.c3*eq.dR2ds*eq.Omega*eq.rhorot)/eq.g - (eq.g12*eq.kappa*eq.Kpar*eq.rhorot*eq.Uthi)/(eq.g*eq.q*eq.R) - (2*eq.c3*eq.g12*eq.Omega*eq.rhorot*imag*ntor)/(eq.g*eq.q)



        

class M_40:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,0,0,1,],[0,0,1,0,],[-3,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (-(eq.c2*eq.dg11du*eq.dR2du*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q**3) + (eq.c2*eq.dg11du*eq.dR2ds*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q**3) - (eq.c2*eq.dR2ds*eq.dR2du*eq.F*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.c2*eq.dR2ds**2*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.c2*eq.dqds*eq.dR2du*eq.F*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.c2*eq.dqds*eq.dR2ds*eq.F*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.B2*eq.c2*eq.dR2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dR2ds**2*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g**3*eq.q) + (eq.c2*eq.dB2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g**3*eq.q) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.q))*imag
        self.var_D00[1] = (eq.dFds*ntor)/eq.g + (eq.dqds*eq.F*ntor)/(eq.g*eq.q) + (eq.dProtds*eq.R2*ntor)/(eq.F*eq.g) - (eq.B2*eq.dqds*eq.R2*ntor)/(eq.F*eq.g*eq.q) - (eq.dR2ds*eq.Prot*eq.R2*eq.U*ntor)/(eq.F*eq.g) + imag*(-((eq.dProtds*eq.dR2du)/(eq.F*eq.g*eq.q)) + (eq.dFds*eq.dR2du*eq.gamma*eq.Prot)/(eq.F**2*eq.g*eq.q) - (eq.dR2ds*eq.dR2du*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q*eq.R2) - (eq.c2*eq.dR2du*eq.drhorotds*eq.Omega**2*eq.R2)/(2.*eq.F*eq.g*eq.q) + (eq.c0*eq.dR2ds*eq.dR2du*eq.Omega**2*eq.rhorot)/(4.*eq.F*eq.g*eq.q) - (eq.c1*eq.dR2ds*eq.dR2du*eq.Omega**2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) - (eq.c2*eq.dOmegads*eq.dR2du*eq.Omega*eq.R2*eq.rhorot)/(eq.F*eq.g*eq.q) + (eq.c1*eq.dFds*eq.dR2du*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F**2*eq.g*eq.q) - (eq.c2*eq.d2R2dsdu*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) + (eq.c2*eq.dg11du*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.F*eq.g*eq.g11*eq.q) - (eq.F*eq.g12*ntor**2)/(eq.g*eq.q*eq.R2) + (eq.c0*eq.g12*eq.Omega**2*eq.R2*eq.rhorot*ntor**2)/(eq.F*eq.g*eq.q))
        self.var_D00[2] = (-(eq.dFds/(eq.g*eq.q)) - (eq.dProtds*eq.R2)/(eq.F*eq.g*eq.q) + (eq.dR2ds*eq.Prot*eq.R2*eq.U)/(eq.F*eq.g*eq.q))*imag - (eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[3] = (-((eq.dR2ds*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q)) - (eq.dProtds*eq.R2)/(eq.F*eq.g*eq.q) + (eq.dFds*eq.gamma*eq.Prot*eq.R2)/(eq.F**2*eq.g*eq.q))*imag

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = (-((eq.dR2du*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q)) - (eq.c1*eq.dR2du*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q))*imag - (eq.F*ntor)/eq.g + (eq.B2*eq.R2*ntor)/(eq.F*eq.g)
        self.var_D01[1] = (-(eq.F/(eq.g*eq.q)) - (eq.gamma*eq.Prot*eq.R2)/(eq.F*eq.g*eq.q))*imag


        

class M_41:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,0,0,1,1,],[0,0,0,1,0,1,],[-3,-1,1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(6,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.c2*eq.dg12du*eq.dR2du*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**3) + (eq.c2*eq.dg12du*eq.dR2ds*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**3) - (eq.c2*eq.dR2du**2*eq.F*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.c2*eq.dR2ds**2*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q) - (eq.c2*eq.dqds*eq.dR2du*eq.F*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.c2*eq.dqds*eq.dR2ds*eq.F*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**2) + (eq.B2*eq.c2*eq.dR2du**2*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dR2ds**2*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.g11*eq.q) + (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.g11*eq.q**2) + (eq.c2*eq.dB2du*eq.dR2du*eq.g11*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g**3*eq.q) - (eq.c2*eq.dB2du*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.q) + (eq.c2*eq.dB2ds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g**3*eq.g11*eq.q) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2**2*eq.rhorot)/(4.*eq.F*eq.g**3*eq.g11*eq.q)
        self.var_D00[1] = -((eq.dProtdu*eq.dR2du)/(eq.F*eq.g*eq.q)) - (eq.dR2du**2*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q*eq.R2) - (eq.c2*eq.dR2du*eq.drhorotdu*eq.Omega**2*eq.R2)/(2.*eq.F*eq.g*eq.q) + (eq.c0*eq.dR2du**2*eq.Omega**2*eq.rhorot)/(4.*eq.F*eq.g*eq.q) - (eq.c1*eq.dR2du**2*eq.Omega**2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) - (eq.c2*eq.d2R2du*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) + (eq.c2*eq.dg12du*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.g11*eq.q) + (eq.c2*eq.dR2ds**2*eq.Omega**2*eq.q*eq.R2*eq.rhorot)/(4.*eq.F*eq.g*eq.g11) + (eq.c2*eq.dqds*eq.dR2ds*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F*eq.g*eq.g11) - (eq.B2*eq.c2*eq.dR2ds**2*eq.Omega**2*eq.q*eq.R2**2*eq.rhorot)/(2.*eq.F**3*eq.g*eq.g11) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.Omega**2*eq.R2**3*eq.rhorot)/(2.*eq.F**3*eq.g*eq.g11) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2**3*eq.rhorot)/(2.*eq.F**4*eq.g*eq.g11) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2**3*eq.rhorot)/(4.*eq.F**3*eq.g*eq.g11) + (eq.F*eq.q*ntor**2)/eq.g - (eq.B2*eq.q*eq.R2*ntor**2)/(eq.F*eq.g) + (eq.c0*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot*ntor**2)/(eq.F*eq.g*eq.g11*eq.q)
        self.var_D00[2] = (eq.c0*eq.g*eq.Omega**2*eq.q*eq.R2**2*eq.rhorot*ntor**2)/(eq.F**3*eq.g11)
        self.var_D00[3] = -((eq.dR2du*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q)) - (eq.c1*eq.dR2du*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.q)
        self.var_D00[4] = -((eq.dR2du*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q)) - (eq.dProtdu*eq.R2)/(eq.F*eq.g*eq.q)
        self.var_D00[5] = -(eq.F/(eq.g*eq.q)) - (eq.gamma*eq.Prot*eq.R2)/(eq.F*eq.g*eq.q)



        

class M_42:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,0,0,1,1,],[0,0,0,1,0,1,],[-3,-1,1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(6,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (-(eq.c2*eq.dqds*eq.dR2du*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.c2*eq.dqds*eq.dR2ds*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**2) + (eq.B2*eq.c2*eq.dR2du**2*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dR2ds**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q) - (eq.c2*eq.dg12du*eq.dR2du*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**3*eq.R2) + (eq.c2*eq.dg12du*eq.dR2ds*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**3*eq.R2) - (eq.c2*eq.dR2du**2*eq.F**2*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q*eq.R2) + (eq.c2*eq.dR2ds**2*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q*eq.R2) + (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**2) + (eq.c2*eq.dB2du*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.c2*eq.dB2du*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.c2*eq.dB2ds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.g11*eq.q))*imag
        self.var_D00[1] = (eq.dR2du*eq.gamma*eq.Prot*ntor)/(eq.g*eq.R2) - (eq.c0*eq.dR2du*eq.Omega**2*eq.rhorot*ntor)/eq.g + (eq.c1*eq.dR2du*eq.Omega**2*eq.rhorot*ntor)/(2.*eq.g) + imag*(-(eq.c2*eq.dR2du*eq.drhorotdu*eq.Omega**2)/(2.*eq.g*eq.q) - (eq.dProtdu*eq.dR2du)/(eq.g*eq.q*eq.R2) - (eq.c2*eq.d2R2du*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q) + (eq.c2*eq.dg12du*eq.dR2ds*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.g11*eq.q) + (eq.c2*eq.dR2ds**2*eq.Omega**2*eq.q*eq.rhorot)/(4.*eq.g*eq.g11) + (eq.c0*eq.dR2du**2*eq.Omega**2*eq.rhorot)/(4.*eq.g*eq.q*eq.R2) + (eq.c2*eq.dqds*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g*eq.g11) - (eq.B2*eq.c2*eq.dR2ds**2*eq.Omega**2*eq.q*eq.R2*eq.rhorot)/(2.*eq.F**2*eq.g*eq.g11) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g*eq.g11) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2**2*eq.rhorot)/(2.*eq.F**3*eq.g*eq.g11) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2**2*eq.rhorot)/(4.*eq.F**2*eq.g*eq.g11) + (eq.c0*eq.g12**2*eq.Omega**2*eq.rhorot*ntor**2)/(eq.g*eq.g11*eq.q))
        self.var_D00[2] = (eq.c0*eq.g*eq.Omega**2*eq.q*eq.R2*eq.rhorot*imag*ntor**2)/(eq.F**2*eq.g11)
        self.var_D00[3] = (-((eq.dR2du*eq.gamma*eq.Prot)/(eq.g*eq.q*eq.R2)) - (eq.c1*eq.dR2du*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q))*imag
        self.var_D00[4] = -((eq.dProtdu*imag)/(eq.g*eq.q)) + (eq.gamma*eq.Prot*ntor)/eq.g
        self.var_D00[5] = -((eq.gamma*eq.Prot*imag)/(eq.g*eq.q))



        

class M_43:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = ((eq.F*eq.g12*eq.kappa*eq.Kpar*eq.rhorot*eq.Uthi)/(eq.B2*eq.g*eq.q*eq.R) - (eq.g12*eq.kappa*eq.Kpar*eq.R2*eq.rhorot*eq.Uthi)/(eq.F*eq.g*eq.q*eq.R))*imag + (2*eq.c3*eq.g12*eq.Omega*eq.R2*eq.rhorot*ntor)/(eq.F*eq.g*eq.q)



        

class M_44:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,],[0,0,],[-1,1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.F*eq.kappa*eq.Kpar*eq.q*eq.R2*eq.rhorot*eq.Uthi)/(eq.B2*eq.g*eq.R)) + (2*eq.kappa*eq.Kpar*eq.q*eq.R2**2*eq.rhorot*eq.Uthi)/(eq.F*eq.g*eq.R) - (eq.B2*eq.kappa*eq.Kpar*eq.q*eq.R2**3*eq.rhorot*eq.Uthi)/(eq.F**3*eq.g*eq.R) - (2*eq.c3*eq.g12**2*eq.Omega*eq.R2*eq.rhorot*imag*ntor)/(eq.F*eq.g*eq.g11*eq.q)
        self.var_D00[1] = (-2*eq.c3*eq.g*eq.Omega*eq.q*eq.R2**2*eq.rhorot*imag*ntor)/(eq.F**3*eq.g11)



        

class M_45:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,],[0,0,],[-1,1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = ((eq.c3*eq.dR2du*eq.Omega*eq.rhorot)/eq.g + (eq.kappa*eq.Kpar*eq.q*eq.R2*eq.rhorot*eq.Uthi)/(eq.g*eq.R) - (eq.B2*eq.kappa*eq.Kpar*eq.q*eq.R2**2*eq.rhorot*eq.Uthi)/(eq.F**2*eq.g*eq.R))*imag + (2*eq.c3*eq.g12**2*eq.Omega*eq.rhorot*ntor)/(eq.g*eq.g11*eq.q)
        self.var_D00[1] = (2*eq.c3*eq.g*eq.Omega*eq.q*eq.R2*eq.rhorot*ntor)/(eq.F**2*eq.g11)



        

class M_50:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,0,0,1,],[0,0,1,0,],[-3,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.c2*eq.dqds*eq.dR2du*eq.F**2*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.c2*eq.dqds*eq.dR2ds*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.B2*eq.c2*eq.dR2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dR2ds**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q) - (eq.c2*eq.dg11du*eq.dR2du*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q**3*eq.R2) + (eq.c2*eq.dg11du*eq.dR2ds*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q**3*eq.R2) - (eq.c2*eq.dR2ds*eq.dR2du*eq.F**2*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q*eq.R2) + (eq.c2*eq.dR2ds**2*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q*eq.R2) + (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.c2*eq.dB2ds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q)
        self.var_D00[1] = -(eq.c2*eq.dR2du*eq.drhorotds*eq.Omega**2)/(2.*eq.g*eq.q) - (eq.c2*eq.dOmegads*eq.dR2du*eq.Omega*eq.rhorot)/(eq.g*eq.q) - (eq.c2*eq.d2R2dsdu*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q) + (eq.c1*eq.dFds*eq.dR2du*eq.Omega**2*eq.rhorot)/(2.*eq.F*eq.g*eq.q) + (eq.c2*eq.dg11du*eq.dR2ds*eq.Omega**2*eq.rhorot)/(4.*eq.g*eq.g11*eq.q) + (eq.c0*eq.dR2ds*eq.dR2du*eq.Omega**2*eq.rhorot)/(4.*eq.g*eq.q*eq.R2) - (eq.c1*eq.dR2ds*eq.dR2du*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q*eq.R2) + (eq.c0*eq.g12*eq.Omega**2*eq.rhorot*ntor**2)/(eq.g*eq.q) + imag*(-((eq.dFds*eq.gamma*eq.Prot*ntor)/(eq.F*eq.g)) + (eq.dR2ds*eq.gamma*eq.Prot*ntor)/(eq.g*eq.R2) - (eq.c0*eq.dR2ds*eq.Omega**2*eq.rhorot*ntor)/eq.g + (eq.dR2ds*eq.Prot*eq.U*ntor)/eq.g)
        self.var_D00[2] = -(eq.dProtds/(eq.g*eq.q)) + (eq.dR2ds*eq.Prot*eq.U)/(eq.g*eq.q)
        self.var_D00[3] = -(eq.dProtds/(eq.g*eq.q)) + (eq.dFds*eq.gamma*eq.Prot)/(eq.F*eq.g*eq.q) - (eq.dR2ds*eq.gamma*eq.Prot)/(eq.g*eq.q*eq.R2)

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -(eq.c1*eq.dR2du*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q) + (eq.gamma*eq.Prot*imag*ntor)/eq.g
        self.var_D01[1] = -((eq.gamma*eq.Prot)/(eq.g*eq.q))


        

class M_51:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,0,0,1,1,],[0,0,0,1,0,1,],[-3,-1,1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(6,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = ((eq.c2*eq.dqds*eq.dR2du*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2) - (eq.c2*eq.dqds*eq.dR2ds*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**2) - (eq.B2*eq.c2*eq.dR2du**2*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dR2ds**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q) + (eq.c2*eq.dg12du*eq.dR2du*eq.F**2*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**3*eq.R2) - (eq.c2*eq.dg12du*eq.dR2ds*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**3*eq.R2) + (eq.c2*eq.dR2du**2*eq.F**2*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q*eq.R2) - (eq.c2*eq.dR2ds**2*eq.F**2*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q*eq.R2) - (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.q**2) + (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**2) - (eq.c2*eq.dB2du*eq.dR2du*eq.g11*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.c2*eq.dB2du*eq.dR2ds*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.c2*eq.dB2ds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g12*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.q) + (eq.c2*eq.dB2ds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q) - (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g**3*eq.g11*eq.q))*imag
        self.var_D00[1] = (eq.dProtdu*ntor)/eq.g + (eq.dR2du*eq.gamma*eq.Prot*ntor)/(eq.g*eq.R2) - (eq.c0*eq.dR2du*eq.Omega**2*eq.rhorot*ntor)/eq.g + imag*((eq.c2*eq.dR2du*eq.drhorotdu*eq.Omega**2)/(2.*eq.g*eq.q) + (eq.c2*eq.d2R2du*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q) - (eq.c2*eq.dg12du*eq.dR2ds*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.g11*eq.q) - (eq.c2*eq.dR2ds**2*eq.Omega**2*eq.q*eq.rhorot)/(4.*eq.g*eq.g11) - (eq.c0*eq.dR2du**2*eq.Omega**2*eq.rhorot)/(4.*eq.g*eq.q*eq.R2) + (eq.c1*eq.dR2du**2*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q*eq.R2) - (eq.c2*eq.dqds*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.g*eq.g11) + (eq.B2*eq.c2*eq.dR2ds**2*eq.Omega**2*eq.q*eq.R2*eq.rhorot)/(2.*eq.F**2*eq.g*eq.g11) + (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.Omega**2*eq.R2**2*eq.rhorot)/(2.*eq.F**2*eq.g*eq.g11) - (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2**2*eq.rhorot)/(2.*eq.F**3*eq.g*eq.g11) + (eq.c2*eq.dB2ds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2**2*eq.rhorot)/(4.*eq.F**2*eq.g*eq.g11) - (eq.c0*eq.g12**2*eq.Omega**2*eq.rhorot*ntor**2)/(eq.g*eq.g11*eq.q))
        self.var_D00[2] = -((eq.c0*eq.g*eq.Omega**2*eq.q*eq.R2*eq.rhorot*imag*ntor**2)/(eq.F**2*eq.g11))
        self.var_D00[3] = (eq.c1*eq.dR2du*eq.Omega**2*eq.rhorot*imag)/(2.*eq.g*eq.q) + (eq.gamma*eq.Prot*ntor)/eq.g
        self.var_D00[4] = (eq.dProtdu/(eq.g*eq.q) + (eq.dR2du*eq.gamma*eq.Prot)/(eq.g*eq.q*eq.R2))*imag
        self.var_D00[5] = (eq.gamma*eq.Prot*imag)/(eq.g*eq.q)



        

class M_52:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,0,0,1,1,],[0,0,0,1,0,1,],[-3,-1,1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(6,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.B2*eq.c2*eq.dqds*eq.dR2du*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**2) + (eq.c2*eq.dB2du*eq.dR2du*eq.F*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) - (eq.B2*eq.c2*eq.dFds*eq.dR2du*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q) - (eq.c2*eq.dB2du*eq.dR2ds*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.c2*eq.dB2ds*eq.dR2du*eq.F*eq.g12*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q) - (eq.c2*eq.dg12du*eq.dR2du*eq.F**3*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**3*eq.R2**2) + (eq.c2*eq.dg12du*eq.dR2ds*eq.F**3*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**3*eq.R2**2) + (eq.c0*eq.dR2du**2*eq.F**3*eq.g11*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.q*eq.R2**2) - (eq.c2*eq.dR2du**2*eq.F**3*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q*eq.R2**2) - (eq.c0*eq.dR2ds*eq.dR2du*eq.F**3*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q*eq.R2**2) + (eq.c2*eq.dR2ds*eq.dR2du*eq.F**3*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q*eq.R2**2) + (eq.c0*eq.dR2ds**2*eq.F**3*eq.g12**2*eq.Omega**2*eq.rhorot)/(4.*eq.g**3*eq.g11*eq.q*eq.R2**2) - (eq.c2*eq.dqds*eq.dR2du*eq.F**3*eq.g12*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q**2*eq.R2) + (eq.c2*eq.dqds*eq.dR2ds*eq.F**3*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q**2*eq.R2) + (eq.B2*eq.c2*eq.dR2du**2*eq.F*eq.g11*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.q*eq.R2) - (eq.B2*eq.c2*eq.dR2ds**2*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot)/(2.*eq.g**3*eq.g11*eq.q*eq.R2)
        self.var_D00[1] = -(eq.c2*eq.dR2du*eq.drhorotdu*eq.F*eq.Omega**2)/(2.*eq.g*eq.q*eq.R2) + (eq.c2*eq.dqds*eq.dR2ds*eq.F*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.g11) - (eq.B2*eq.c2*eq.dR2ds**2*eq.Omega**2*eq.q*eq.rhorot)/(2.*eq.F*eq.g*eq.g11) + (eq.c0*eq.dR2du**2*eq.F*eq.Omega**2*eq.rhorot)/(4.*eq.g*eq.q*eq.R2**2) - (eq.c2*eq.d2R2du*eq.F*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q*eq.R2) + (eq.c2*eq.dg12du*eq.dR2ds*eq.F*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.g11*eq.q*eq.R2) + (eq.c0*eq.dR2ds**2*eq.F*eq.Omega**2*eq.q*eq.rhorot)/(4.*eq.g*eq.g11*eq.R2) - (eq.B2*eq.c2*eq.dqds*eq.dR2ds*eq.Omega**2*eq.R2*eq.rhorot)/(2.*eq.F*eq.g*eq.g11) + (eq.B2*eq.c2*eq.dFds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2*eq.rhorot)/(2.*eq.F**2*eq.g*eq.g11) - (eq.c2*eq.dB2ds*eq.dR2ds*eq.Omega**2*eq.q*eq.R2*eq.rhorot)/(4.*eq.F*eq.g*eq.g11) - (eq.F*eq.gamma*eq.Prot*eq.q*ntor**2)/(eq.g*eq.R2) + (eq.c0*eq.F*eq.Omega**2*eq.q*eq.rhorot*ntor**2)/eq.g + (eq.c0*eq.F*eq.g12**2*eq.Omega**2*eq.rhorot*ntor**2)/(eq.g*eq.g11*eq.q*eq.R2) + imag*((eq.dProtdu*eq.F*ntor)/(eq.g*eq.R2) - (eq.c1*eq.dR2du*eq.F*eq.Omega**2*eq.rhorot*ntor)/(2.*eq.g*eq.R2))
        self.var_D00[2] = (eq.c0*eq.g*eq.Omega**2*eq.q*eq.rhorot*ntor**2)/(eq.F*eq.g11)
        self.var_D00[3] = -(eq.c1*eq.dR2du*eq.F*eq.Omega**2*eq.rhorot)/(2.*eq.g*eq.q*eq.R2) + (eq.F*eq.gamma*eq.Prot*imag*ntor)/(eq.g*eq.R2)
        self.var_D00[4] = -((eq.dProtdu*eq.F)/(eq.g*eq.q*eq.R2)) - (eq.F*eq.gamma*eq.Prot*imag*ntor)/(eq.g*eq.R2)
        self.var_D00[5] = -((eq.F*eq.gamma*eq.Prot)/(eq.g*eq.q*eq.R2))



        

class M_53:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,],[-1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.c3*eq.dR2ds*eq.Omega*eq.rhorot)/eq.g) - (eq.g12*eq.kappa*eq.Kpar*eq.rhorot*eq.Uthi)/(eq.g*eq.q*eq.R) - (2*eq.c3*eq.g12*eq.Omega*eq.rhorot*imag*ntor)/(eq.g*eq.q)



        

class M_54:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,],[0,0,],[-1,1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = ((eq.c3*eq.dR2du*eq.Omega*eq.rhorot)/eq.g - (eq.kappa*eq.Kpar*eq.q*eq.R2*eq.rhorot*eq.Uthi)/(eq.g*eq.R) + (eq.B2*eq.kappa*eq.Kpar*eq.q*eq.R2**2*eq.rhorot*eq.Uthi)/(eq.F**2*eq.g*eq.R))*imag - (2*eq.c3*eq.g12**2*eq.Omega*eq.rhorot*ntor)/(eq.g*eq.g11*eq.q)
        self.var_D00[1] = (-2*eq.c3*eq.g*eq.Omega*eq.q*eq.R2*eq.rhorot*ntor)/(eq.F**2*eq.g11)



        

class M_55:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,],[0,0,],[-1,1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.B2*eq.kappa*eq.Kpar*eq.q*eq.R2*eq.rhorot*eq.Uthi)/(eq.F*eq.g*eq.R)) + imag*((-2*eq.c3*eq.F*eq.Omega*eq.q*eq.rhorot*ntor)/eq.g - (2*eq.c3*eq.F*eq.g12**2*eq.Omega*eq.rhorot*ntor)/(eq.g*eq.g11*eq.q*eq.R2))
        self.var_D00[1] = (-2*eq.c3*eq.g*eq.Omega*eq.q*eq.rhorot*imag*ntor)/(eq.F*eq.g11)



        




# Boundary Conditions 
#==============================================================#
class BC_0:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Dirichlet']*grid.Mtot)


class BC_1:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)

        idx = np.where(abs(grid.m) == 1)[0]
        self.Axis[idx] = ['Neumann']

class BC_2:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)

        idx = np.where(abs(grid.m) == 0)[0]
        self.Axis[idx] = ['Neumann']

class BC_3:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Dirichlet']*grid.Mtot)


class BC_4:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)

        idx = np.where(abs(grid.m) == 1)[0]
        self.Axis[idx] = ['Neumann']

class BC_5:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)

        idx = np.where(abs(grid.m) == 0)[0]
        self.Axis[idx] = ['Neumann']


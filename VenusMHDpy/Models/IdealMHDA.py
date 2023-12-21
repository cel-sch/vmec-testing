import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,True,True]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dqds**2*eq.F)/(eq.g*eq.q**3) + (eq.dFds*eq.dqds)/(eq.g*eq.q**2) - (eq.dPds*eq.dR2ds)/(eq.F*eq.g*eq.q) + (2*eq.dFds*eq.dR2ds*eq.gamma*eq.P)/(eq.F**2*eq.g*eq.q) - (eq.dR2ds*eq.dR2du*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**3) + (eq.dg12du*eq.dR2ds*eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2) + (eq.dFds*eq.dR2du*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2) - (eq.dB2du*eq.dR2ds*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q**3*eq.R2**2) - (eq.dFds*eq.dg12du*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2) + (eq.dB2du*eq.dFds*eq.g12*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q**3*eq.R2) - (eq.dR2ds**2*eq.gamma*eq.P)/(eq.F*eq.g*eq.q*eq.R2) - (eq.B2*eq.dqds**2*eq.R2)/(eq.F*eq.g*eq.q**3) + (eq.dPds*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) + (eq.dFds*eq.dPds*eq.R2)/(eq.F**2*eq.g*eq.q) - (eq.dFds**2*eq.gamma*eq.P*eq.R2)/(eq.F**3*eq.g*eq.q) - (eq.F*eq.g11*ntor**2)/(eq.g*eq.q*eq.R2) + imag*((eq.dR2ds*eq.F*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2**2) - (eq.dFds*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2))
        self.var_D00[1] = (eq.dR2ds*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2) - (eq.dqds*eq.F*eq.g12)/(eq.g*eq.q**4*eq.R2) - (eq.dFds*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2) + (eq.F*eq.g11*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[2] = -((eq.dqds*eq.F*eq.g12)/(eq.g*eq.q**4*eq.R2)) - (eq.F*eq.g11*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[3] = -((eq.F*eq.g11)/(eq.g*eq.q**3*eq.R2))

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.dqds*eq.F)/(eq.g*eq.q**2)) - (eq.dR2ds*eq.gamma*eq.P)/(eq.F*eq.g*eq.q) + (eq.B2*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) - (eq.dPds*eq.R2)/(eq.F*eq.g*eq.q) + (eq.dFds*eq.gamma*eq.P*eq.R2)/(eq.F**2*eq.g*eq.q) - (eq.F*eq.g12*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D01[1] = (eq.F*eq.g12)/(eq.g*eq.q**3*eq.R2)

        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.dqds*eq.F)/(eq.g*eq.q**2)) - (eq.dR2ds*eq.gamma*eq.P)/(eq.F*eq.g*eq.q) - (eq.dR2du*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2) + (eq.dg12du*eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2) - (eq.dB2du*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q**3*eq.R2) + (eq.B2*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) - (eq.dPds*eq.R2)/(eq.F*eq.g*eq.q) + (eq.dFds*eq.gamma*eq.P*eq.R2)/(eq.F**2*eq.g*eq.q) + imag*((eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2) + (eq.F*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2))
        self.var_D10[1] = (eq.F*eq.g12)/(eq.g*eq.q**3*eq.R2) + (eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2)

        self.DerPol_DerBs11 = [[0,],[0,],[-1,]]

        self.var_D11 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D11[0] = -((eq.B2*eq.R2)/(eq.F*eq.g*eq.q)) - (eq.gamma*eq.P*eq.R2)/(eq.F*eq.g*eq.q)
        

class M_01:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,0,1,],[0,1,0,],[-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(3,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dFds*ntor)/eq.g + (eq.dFds*eq.gamma*eq.P*ntor)/(eq.B2*eq.g) + (eq.dR2ds*eq.gamma*eq.P*ntor)/(eq.F*eq.g) + (eq.dqds*eq.F*ntor)/(eq.g*eq.q) - (eq.dR2ds*eq.F*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.R2) + (eq.dPds*eq.R2*ntor)/(eq.F*eq.g) - (eq.dFds*eq.gamma*eq.P*eq.R2*ntor)/(eq.F**2*eq.g) - (eq.B2*eq.dqds*eq.R2*ntor)/(eq.F*eq.g*eq.q) + imag*((eq.dB2du*eq.dFds*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q) - (eq.dB2du*eq.dR2ds*eq.F*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q*eq.R2) + (eq.F*eq.g12*ntor**2)/(eq.g*eq.q*eq.R2))
        self.var_D00[1] = (-(eq.dFds/(eq.g*eq.q)) - (eq.dFds*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q) + (eq.dR2ds*eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q*eq.R2))*imag
        self.var_D00[2] = -((eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2))


        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.dB2du*eq.F*eq.gamma*eq.P*imag)/(eq.B2**2*eq.g*eq.q)) - (eq.F*ntor)/eq.g - (eq.F*eq.gamma*eq.P*ntor)/(eq.B2*eq.g) + (eq.B2*eq.R2*ntor)/(eq.F*eq.g) + (eq.gamma*eq.P*eq.R2*ntor)/(eq.F*eq.g)
        self.var_D10[1] = (eq.F/(eq.g*eq.q) + (eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q))*imag

        

class M_02:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dFds*eq.gamma*eq.P*ntor)/(eq.F*eq.g) - (eq.dR2ds*eq.gamma*eq.P*ntor)/(eq.g*eq.R2)
        self.var_D00[1] = (-((eq.dFds*eq.gamma*eq.P)/(eq.F*eq.g*eq.q)) + (eq.dR2ds*eq.gamma*eq.P)/(eq.g*eq.q*eq.R2))*imag


        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.gamma*eq.P*ntor)/eq.g)
        self.var_D10[1] = (eq.gamma*eq.P*imag)/(eq.g*eq.q)

        

class M_10:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dFds*ntor)/eq.g + (eq.dqds*eq.F*ntor)/(eq.g*eq.q) - (eq.dR2du*eq.F*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2**2) + (eq.dPds*eq.R2*ntor)/(eq.F*eq.g) - (eq.B2*eq.dqds*eq.R2*ntor)/(eq.F*eq.g*eq.q) + imag*((eq.dFds*eq.dR2du*eq.gamma*eq.P)/(eq.F**2*eq.g*eq.q) - (eq.dR2du**2*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**3) + (eq.dg12du*eq.dR2du*eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2) - (eq.dB2du*eq.dR2du*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q**3*eq.R2**2) - (eq.dR2ds*eq.dR2du*eq.gamma*eq.P)/(eq.F*eq.g*eq.q*eq.R2) - (eq.F*eq.g12*ntor**2)/(eq.g*eq.q*eq.R2))
        self.var_D00[1] = (-(eq.dFds/(eq.g*eq.q)) + (eq.dR2du*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2))*imag - (eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[2] = (-((eq.dR2ds*eq.gamma*eq.P)/(eq.F*eq.g*eq.q)) - (eq.dR2du*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2) + (eq.dg12du*eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2) - (eq.dB2du*eq.F*eq.g12*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q**3*eq.R2) + (eq.dFds*eq.gamma*eq.P*eq.R2)/(eq.F**2*eq.g*eq.q))*imag - (eq.F*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2)
        self.var_D00[3] = (eq.F*eq.g12*eq.gamma*eq.P*imag)/(eq.B2*eq.g*eq.q**3*eq.R2)

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.dR2du*eq.gamma*eq.P*imag)/(eq.F*eq.g*eq.q)) - (eq.F*ntor)/eq.g + (eq.B2*eq.R2*ntor)/(eq.F*eq.g)
        self.var_D01[1] = (-(eq.F/(eq.g*eq.q)) - (eq.gamma*eq.P*eq.R2)/(eq.F*eq.g*eq.q))*imag


        

class M_11:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dB2du*eq.dR2du*eq.F*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q*eq.R2) + (eq.F*eq.q*ntor**2)/eq.g - (eq.B2*eq.q*eq.R2*ntor**2)/(eq.F*eq.g) + imag*((eq.dR2du*eq.gamma*eq.P*ntor)/(eq.F*eq.g) - (eq.dR2du*eq.F*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.R2))
        self.var_D00[1] = -((eq.dR2du*eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q*eq.R2))
        self.var_D00[2] = (eq.dB2du*eq.F*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q) + imag*(-((eq.F*eq.gamma*eq.P*ntor)/(eq.B2*eq.g)) + (eq.gamma*eq.P*eq.R2*ntor)/(eq.F*eq.g))
        self.var_D00[3] = -(eq.F/(eq.g*eq.q)) - (eq.F*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q)



        

class M_12:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.dR2du*eq.gamma*eq.P*imag*ntor)/(eq.g*eq.R2))
        self.var_D00[1] = -((eq.dR2du*eq.gamma*eq.P)/(eq.g*eq.q*eq.R2))
        self.var_D00[2] = -((eq.gamma*eq.P*imag*ntor)/eq.g)
        self.var_D00[3] = -((eq.gamma*eq.P)/(eq.g*eq.q))



        

class M_20:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dFds*eq.gamma*eq.P*ntor)/(eq.F*eq.g) - (eq.dR2du*eq.F**2*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2**3) + (eq.dg12du*eq.F**2*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2**2) - (eq.dB2du*eq.F**2*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2**2*eq.g*eq.q**2*eq.R2**2) - (eq.dR2ds*eq.gamma*eq.P*ntor)/(eq.g*eq.R2) + (eq.F**2*eq.g12*eq.gamma*eq.P*imag*ntor**2)/(eq.B2*eq.g*eq.q*eq.R2**2)
        self.var_D00[1] = (eq.F**2*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2**2)
        self.var_D00[2] = ((eq.dFds*eq.gamma*eq.P)/(eq.F*eq.g*eq.q) - (eq.dR2du*eq.F**2*eq.g12*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**3) + (eq.dg12du*eq.F**2*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q**3*eq.R2**2) - (eq.dB2du*eq.F**2*eq.g12*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q**3*eq.R2**2) - (eq.dR2ds*eq.gamma*eq.P)/(eq.g*eq.q*eq.R2))*imag - (eq.F**2*eq.g12*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.q**2*eq.R2**2)
        self.var_D00[3] = (eq.F**2*eq.g12*eq.gamma*eq.P*imag)/(eq.B2*eq.g*eq.q**3*eq.R2**2)

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.gamma*eq.P*ntor)/eq.g)
        self.var_D01[1] = -((eq.gamma*eq.P*imag)/(eq.g*eq.q))


        

class M_21:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.dB2du*eq.F**2*eq.gamma*eq.P*imag*ntor)/(eq.B2**2*eq.g*eq.R2)) + (eq.gamma*eq.P*eq.q*ntor**2)/eq.g - (eq.F**2*eq.gamma*eq.P*eq.q*ntor**2)/(eq.B2*eq.g*eq.R2)
        self.var_D00[1] = (eq.F**2*eq.gamma*eq.P*imag*ntor)/(eq.B2*eq.g*eq.R2)
        self.var_D00[2] = (eq.dB2du*eq.F**2*eq.gamma*eq.P)/(eq.B2**2*eq.g*eq.q*eq.R2) + imag*((eq.gamma*eq.P*ntor)/eq.g - (eq.F**2*eq.gamma*eq.P*ntor)/(eq.B2*eq.g*eq.R2))
        self.var_D00[3] = -((eq.F**2*eq.gamma*eq.P)/(eq.B2*eq.g*eq.q*eq.R2))



        

class M_22:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.F*eq.gamma*eq.P*eq.q*ntor**2)/(eq.g*eq.R2))
        self.var_D00[1] = (eq.F*eq.gamma*eq.P*imag*ntor)/(eq.g*eq.R2)
        self.var_D00[2] = -((eq.F*eq.gamma*eq.P*imag*ntor)/(eq.g*eq.R2))
        self.var_D00[3] = -((eq.F*eq.gamma*eq.P)/(eq.g*eq.q*eq.R2))


# Boundary conditions
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
        


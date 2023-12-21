import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,True,True]

        self.DerPol_DerBs00 = [[0,0,1,1,],[0,1,0,1,],[-1,-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(4,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dqds**2*eq.F)/(eq.g*eq.q**3) + (eq.dFds*eq.dqds)/(eq.g*eq.q**2) - (eq.dPds*eq.dR2ds)/(eq.F*eq.g*eq.q) - (eq.B2*eq.dqds**2*eq.R2)/(eq.F*eq.g*eq.q**3) + (eq.dPds*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) + (eq.dFds*eq.dPds*eq.R2)/(eq.F**2*eq.g*eq.q) - (eq.F*eq.g11*ntor**2)/(eq.g*eq.q*eq.R2)
        self.var_D00[1] = -((eq.dqds*eq.F*eq.g12)/(eq.g*eq.q**4*eq.R2)) + (eq.F*eq.g11*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[2] = -((eq.dqds*eq.F*eq.g12)/(eq.g*eq.q**4*eq.R2)) - (eq.F*eq.g11*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D00[3] = -((eq.F*eq.g11)/(eq.g*eq.q**3*eq.R2))

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.dqds*eq.F)/(eq.g*eq.q**2)) + (eq.B2*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) - (eq.dPds*eq.R2)/(eq.F*eq.g*eq.q) - (eq.F*eq.g12*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D01[1] = (eq.F*eq.g12)/(eq.g*eq.q**3*eq.R2)

        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.dqds*eq.F)/(eq.g*eq.q**2)) + (eq.B2*eq.dqds*eq.R2)/(eq.F*eq.g*eq.q**2) - (eq.dPds*eq.R2)/(eq.F*eq.g*eq.q) + (eq.F*eq.g12*imag*ntor)/(eq.g*eq.q**2*eq.R2)
        self.var_D10[1] = (eq.F*eq.g12)/(eq.g*eq.q**3*eq.R2)

        self.DerPol_DerBs11 = [[0,],[0,],[-1,]]

        self.var_D11 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D11[0] = -((eq.B2*eq.R2)/(eq.F*eq.g*eq.q))
        

class M_01:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,0,1,],[0,1,0,],[-1,-1,-1,]]

        self.var_D00 = np.zeros(shape=(3,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dFds*ntor)/eq.g + (eq.dqds*eq.F*ntor)/(eq.g*eq.q) + (eq.dPds*eq.R2*ntor)/(eq.F*eq.g) - (eq.B2*eq.dqds*eq.R2*ntor)/(eq.F*eq.g*eq.q) + (eq.F*eq.g12*imag*ntor**2)/(eq.g*eq.q*eq.R2)
        self.var_D00[1] = -((eq.dFds*imag)/(eq.g*eq.q))
        self.var_D00[2] = -((eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2))


        self.DerPol_DerBs10 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.F*ntor)/eq.g) + (eq.B2*eq.R2*ntor)/(eq.F*eq.g)
        self.var_D10[1] = (eq.F*imag)/(eq.g*eq.q)

        

class M_10:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,0,],[0,1,],[-1,-1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dFds*ntor)/eq.g + (eq.dqds*eq.F*ntor)/(eq.g*eq.q) + (eq.dPds*eq.R2*ntor)/(eq.F*eq.g) - (eq.B2*eq.dqds*eq.R2*ntor)/(eq.F*eq.g*eq.q) - (eq.F*eq.g12*imag*ntor**2)/(eq.g*eq.q*eq.R2)
        self.var_D00[1] = -((eq.dFds*imag)/(eq.g*eq.q)) - (eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2)

        self.DerPol_DerBs01 = [[0,1,],[0,0,],[-1,-1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.F*ntor)/eq.g) + (eq.B2*eq.R2*ntor)/(eq.F*eq.g)
        self.var_D01[1] = -((eq.F*imag)/(eq.g*eq.q))


        

class M_11:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,1,],[-1,-1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.F*eq.q*ntor**2)/eq.g - (eq.B2*eq.q*eq.R2*ntor**2)/(eq.F*eq.g)
        self.var_D00[1] = -(eq.F/(eq.g*eq.q))



        




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


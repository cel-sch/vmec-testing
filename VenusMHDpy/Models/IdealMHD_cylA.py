import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,True,True]

        self.DerPol_DerBs00 = [[0,0,0,],[0,1,2,]]

        self.var_D00 = np.zeros(shape=(3,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz**2*eq.dqds*eq.S)/eq.q**3 - (eq.Bz*eq.dBzds*eq.S)/eq.q**2 + (eq.Bz**2*k**2)/2.
        self.var_D00[1] = -((eq.Bz**2*imag*k)/eq.q)
        self.var_D00[2] = -eq.Bz**2/(2.*eq.q**2)

        self.DerPol_DerBs01 = [[0,],[0,]]

        self.var_D01 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.Bz**2*eq.S)/eq.q**2)

        self.DerPol_DerBs10 = [[0,],[0,]]

        self.var_D10 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.Bz**2*eq.S)/eq.q**2)

        self.DerPol_DerBs11 = [[0,],[0,]]

        self.var_D11 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D11[0] = eq.Bz**2/2. + (eq.gamma*eq.P)/2. + (eq.Bz**2*eq.S**2)/(2.*eq.q**2)
        

class M_01:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz**2*eq.S*k)/eq.q


        self.DerPol_DerBs10 = [[0,0,],[0,1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -(eq.Bz**2*eq.S**2*k)/(2.*eq.q)
        self.var_D10[1] = (-eq.Bz**2/2. - (eq.gamma*eq.P)/2.)*imag

        

class M_02:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,True,False]



        self.DerPol_DerBs10 = [[0,0,],[0,1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = (eq.Bz*eq.gamma*eq.P*k)/2.
        self.var_D10[1] = -(eq.Bz*eq.gamma*eq.P*imag)/(2.*eq.q)

        

class M_10:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz**2*eq.S*k)/eq.q

        self.DerPol_DerBs01 = [[0,0,],[0,1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -(eq.Bz**2*eq.S**2*k)/(2.*eq.q)
        self.var_D01[1] = (-eq.Bz**2/2. - (eq.gamma*eq.P)/2.)*imag


        

class M_11:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,],[0,2,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz**2*eq.S**2*k**2)/2.
        self.var_D00[1] = -eq.Bz**2/2. - (eq.gamma*eq.P)/2.



        

class M_12:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,],[1,2,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.Bz*eq.gamma*eq.P*imag*k)/2.
        self.var_D00[1] = -(eq.Bz*eq.gamma*eq.P)/(2.*eq.q)



        

class M_20:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,True,False,False]


        self.DerPol_DerBs01 = [[0,0,],[0,1,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = (eq.Bz*eq.gamma*eq.P*k)/2.
        self.var_D01[1] = -(eq.Bz*eq.gamma*eq.P*imag)/(2.*eq.q)


        

class M_21:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,],[1,2,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -(eq.Bz*eq.gamma*eq.P*imag*k)/2.
        self.var_D00[1] = -(eq.Bz*eq.gamma*eq.P)/(2.*eq.q)



        

class M_22:
    def __init__(self,eq,k):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,0,0,],[0,1,2,]]

        self.var_D00 = np.zeros(shape=(3,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.Bz**2*eq.gamma*eq.P*k**2)/2.
        self.var_D00[1] = -((eq.Bz**2*eq.gamma*eq.P*imag*k)/eq.q)
        self.var_D00[2] = -(eq.Bz**2*eq.gamma*eq.P)/(2.*eq.q**2)


# Boundary conditions
#==============================================================#


class BC_0:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Dirichlet']*grid.Mtot)
        
class BC_1:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)
        
        idx = np.where(abs(grid.m) == 1)[0]
        self.Axis[idx] = ['Neumann']
        
class BC_2:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)
        
        idx = np.where(abs(grid.m) == 0)[0]
        self.Axis[idx] = ['Neumann']
        


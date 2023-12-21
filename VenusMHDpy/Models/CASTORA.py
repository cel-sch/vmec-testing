import numpy as np 

imag = 1.0j 

class M_00:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_01:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.drhods*eq.R2)/eq.g) - (eq.dR2ds*eq.rho)/eq.g

        self.DerPol_DerBs01 = [[0,],[0,]]

        self.var_D01 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.R2*eq.rho)/eq.g)


        

class M_02:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[1,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.R2*eq.rho*imag)/eq.g)



        

class M_03:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.q*eq.R2*eq.rho*ntor)/eq.g)
        self.var_D00[1] = -((eq.R2*eq.rho*imag)/eq.g)



        

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
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dR2ds*eq.T)/eq.g


        self.DerPol_DerBs10 = [[0,],[0,]]

        self.var_D10 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = (eq.R2*eq.T)/eq.g

        

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
        self.DerBs = [False,False,False,False]




        

class M_14:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dR2ds*eq.rho)/eq.g


        self.DerPol_DerBs10 = [[0,],[0,]]

        self.var_D10 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = (eq.R2*eq.rho)/eq.g

        

class M_15:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,False]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dqds*eq.F**2*ntor)/(eq.g*eq.q) + (eq.dPds*eq.R2*ntor)/eq.g + (eq.B2*eq.dFds*eq.R2*ntor)/(eq.F*eq.g) - (eq.B2*eq.dqds*eq.R2*ntor)/(eq.g*eq.q) + (eq.F**2*eq.g12*imag*ntor**2)/(eq.g*eq.q*eq.R2)
        self.var_D00[1] = -((eq.F**2*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2))


        self.DerPol_DerBs10 = [[0,0,],[0,1,]]

        self.var_D10 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.F**2*ntor)/eq.g) + (eq.B2*eq.R2*ntor)/eq.g
        self.var_D10[1] = (eq.F**2*imag)/(eq.g*eq.q)

        

class M_16:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,True,True]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.F**2*eq.g11*ntor**2)/(eq.g*eq.q*eq.R2) + imag*((eq.dqds*eq.F**2*eq.g12*ntor)/(eq.g*eq.q**3*eq.R2) - (eq.dFds*eq.F*eq.g12*ntor)/(eq.g*eq.q**2*eq.R2))
        self.var_D00[1] = (eq.F**2*eq.g11*imag*ntor)/(eq.g*eq.q**2*eq.R2)


        self.DerPol_DerBs10 = [[0,],[0,]]

        self.var_D10 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = -((eq.F**2*eq.g12*imag*ntor)/(eq.g*eq.q**2*eq.R2))

        self.DerPol_DerBs11 = [[0,],[0,]]

        self.var_D11 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D11[0] = eq.F**2/(eq.g*eq.q)
        

class M_17:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,True,True]

        self.DerPol_DerBs00 = [[0,1,],[1,1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.dqds*eq.F**2*eq.g12)/(eq.g*eq.q**3*eq.R2)) + (eq.dFds*eq.F*eq.g12)/(eq.g*eq.q**2*eq.R2) + (eq.F**2*eq.g11*imag*ntor)/(eq.g*eq.q*eq.R2)
        self.var_D00[1] = -((eq.F**2*eq.g11)/(eq.g*eq.q**2*eq.R2))

        self.DerPol_DerBs01 = [[0,1,],[0,0,]]

        self.var_D01 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.dqds*eq.F**2)/(eq.g*eq.q)) - (eq.dPds*eq.R2)/eq.g - (eq.B2*eq.dFds*eq.R2)/(eq.F*eq.g) + (eq.B2*eq.dqds*eq.R2)/(eq.g*eq.q) - (eq.F**2*eq.g12*imag*ntor)/(eq.g*eq.q*eq.R2)
        self.var_D01[1] = (eq.F**2*eq.g12)/(eq.g*eq.q**2*eq.R2)

        self.DerPol_DerBs10 = [[0,],[1,]]

        self.var_D10 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D10[0] = (eq.F**2*eq.g12)/(eq.g*eq.q**2*eq.R2)

        self.DerPol_DerBs11 = [[0,],[0,]]

        self.var_D11 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D11[0] = eq.F**2/eq.g - (eq.B2*eq.R2)/eq.g
        

class M_20:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dR2du*eq.T*imag)/eq.g
        self.var_D00[1] = (eq.R2*eq.T*imag)/eq.g



        

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
        self.DerBs = [False,False,False,False]




        

class M_24:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dR2du*eq.rho*imag)/eq.g
        self.var_D00[1] = (eq.R2*eq.rho*imag)/eq.g



        

class M_25:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,1,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.F**2*eq.q*ntor**2)/eq.g - (eq.B2*eq.q*eq.R2*ntor**2)/eq.g
        self.var_D00[1] = -(eq.F**2/(eq.g*eq.q))



        

class M_26:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.dFds*eq.F*ntor)/eq.g) - (eq.dPds*eq.R2*ntor)/eq.g + (eq.F**2*eq.g12*imag*ntor**2)/(eq.g*eq.q*eq.R2)

        self.DerPol_DerBs01 = [[1,],[0,]]

        self.var_D01 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = (eq.F**2*imag)/(eq.g*eq.q)


        

class M_27:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,],[1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (-((eq.dFds*eq.F)/eq.g) - (eq.dPds*eq.R2)/eq.g)*imag - (eq.F**2*eq.g12*ntor)/(eq.g*eq.q*eq.R2)

        self.DerPol_DerBs01 = [[0,],[0,]]

        self.var_D01 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.F**2*eq.q*ntor)/eq.g) + (eq.B2*eq.q*eq.R2*ntor)/eq.g


        

class M_30:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dR2du*eq.T*imag)/eq.g + (eq.q*eq.R2*eq.T*ntor)/eq.g
        self.var_D00[1] = (eq.R2*eq.T*imag)/eq.g



        

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
        self.DerBs = [False,False,False,False]




        

class M_34:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.dR2du*eq.rho*imag)/eq.g + (eq.q*eq.R2*eq.rho*ntor)/eq.g
        self.var_D00[1] = (eq.R2*eq.rho*imag)/eq.g



        

class M_35:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_36:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.dPds*eq.R2*ntor)/eq.g)



        

class M_37:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[1,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.dPds*eq.R2*imag)/eq.g)



        

class M_40:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

class M_41:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,True,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.dTds*eq.R2*eq.rho)/(eq.g*(-1 + eq.gamma))) + (eq.dR2ds*eq.rho*eq.T)/(eq.g*(-1 + eq.gamma)) - (eq.dR2ds*eq.gamma*eq.rho*eq.T)/(eq.g*(-1 + eq.gamma))

        self.DerPol_DerBs01 = [[0,],[0,]]

        self.var_D01 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D01[0] = -((eq.R2*eq.rho*eq.T)/eq.g)


        

class M_42:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[1,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.R2*eq.rho*eq.T*imag)/eq.g)



        

class M_43:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,1,],[0,0,]]

        self.var_D00 = np.zeros(shape=(2,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.q*eq.R2*eq.rho*eq.T*ntor)/eq.g)
        self.var_D00[1] = -((eq.R2*eq.rho*eq.T*imag)/eq.g)



        

class M_44:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [False,False,False,False]




        

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
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.F**2*eq.g12*imag)/(eq.g*eq.q)



        

class M_52:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.F**2*eq.g12**2)/(eq.g*eq.g11*eq.q) + (eq.g*eq.h**2*eq.q*eq.R2)/eq.g11



        

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
        self.DerBs = [False,False,False,False]




        

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
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = -((eq.F**2*eq.g11)/(eq.g*eq.q))



        

class M_62:
    def __init__(self,eq,ntor):

        #Convention is: 00,01,10,11.
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = (eq.F**2*eq.g12*imag)/(eq.g*eq.q)



        

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
        self.DerBs = [False,False,False,False]




        

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
        self.DerBs = [True,False,False,False]

        self.DerPol_DerBs00 = [[0,],[0,]]

        self.var_D00 = np.zeros(shape=(1,eq.R.shape[0], eq.R.shape[1]), dtype=complex)
        self.var_D00[0] = eq.g*eq.h**2



        

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
        self.DerBs = [False,False,False,False]


# Boundary conditions
#==============================================================#


class BC_0:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)
        
        idx = np.where(abs(grid.m) == 0)[0]
        self.Axis[idx] = ['Neumann']
        
class BC_1:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Dirichlet']*grid.Mtot)
        
class BC_2:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)
        
        idx = np.where(abs(grid.m) == 1)[0]
        self.Axis[idx] = ['Neumann']

class BC_3:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)
        
        idx = np.where(abs(grid.m) == 0)[0]
        self.Axis[idx] = ['Neumann']
        
class BC_4:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)
        
        idx = np.where(abs(grid.m) == 0)[0]
        self.Axis[idx] = ['Neumann']

class BC_5:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Natural']*grid.Mtot)
        self.Edge = np.asarray(['Natural']*grid.Mtot)
        
        idx = np.where(abs(grid.m) == 1)[0]
        self.Axis[idx] = ['Neumann']

class BC_6:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Dirichlet']*grid.Mtot)
        
class BC_7:
    def __init__(self,grid,eq,ntor):

        self.Axis = np.asarray(['Dirichlet']*grid.Mtot)
        self.Edge = np.asarray(['Dirichlet']*grid.Mtot)


        


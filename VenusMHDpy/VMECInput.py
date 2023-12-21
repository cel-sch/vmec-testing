import numpy as np


#*****************************************************************************
#**                        READ VMEC INPUT FILE                             **  
#**                                                                         **
#** This script reads a VMEC input file and stores the variables in a class **
#**                                                                         **
#** AUTHOR: Guillermo Bustos Ramirez                                        ** 
#** LAST CHANGE: 18/09/2019                                                 **
#*****************************************************************************

class RawData:
    def __init__(self, File):
        F = open(File,'r')
        lines = F.readlines()
        F.close()

        str = 'name'
        variable = ''
        value = ''
        
        for i in range(1,len(lines),1):
            for j in lines[i]:
                
                if j=='!':
                    break
                
                if str=='name':
                    if j!='=':
                        variable += j
                    if j=='=':
                        str = 'value'
                
                if str=='value':
                    if j!=',':
                        value += j
                    if j==',':
                        str = 'name'
            
                        varname = variable.upper().strip()
                        vars(self)[varname] = value.split('=')[1].strip().split()
                        
                        #Exceptions for variables with the same name and parenthesis RBC(0,0), EXTCUR(0), etc.
                        exeptions_Boundary = ['RBC', 'RBS', 'ZBC', 'ZBS']
                        
                        for s in exeptions_Boundary:
                            if s in varname:
                                n = varname.split('(')[1].split(',')[0]
                                m = varname.split(',')[1].split(')')[0]
                                varname = varname.split('(')[0]
                                
                                if hasattr(self, varname):
                                    vars(self)[varname].append(float(value.split('=')[1].strip().split()[0]))
                                else:
                                    vars(self)[varname] = [float(value.split('=')[1].strip().split()[0])]
                                
                                if hasattr(self, varname+'n'): 
                                    vars(self)[varname+'n'].append(int(n))
                                else:
                                    vars(self)[varname+'n'] = [int(n)]
                                    
                                if hasattr(self, varname+'m'):
                                    vars(self)[varname+'m'].append(int(m))
                                else:
                                    vars(self)[varname+'m'] = [int(m)]
                                    
                        if 'EXTCUR' in varname:
                            c = varname.split('(')[1].split(')')[0]
                            varname = varname.split('(')[0]
                            
                            if hasattr(self, varname):
                                vars(self)[varname].append(float(value.split('=')[1].strip().split()[0]))
                            else:
                                vars(self)[varname] = [float(value.split('=')[1].strip().split()[0])]
                            
                            if hasattr(self, varname+'c'): 
                                vars(self)[varname+'c'].append(int(c))
                            else:
                                vars(self)[varname+'c'] = [int(c)]
            
                        variable = ''
                        value = ''
                        
def GetVar(RawClass, varname, vartype):
    
    try:
        if vartype.lower() == 'str':
            return vars(RawClass)[varname][0]
        if vartype.lower() == 'int':
            return int(vars(RawClass)[varname][0])
        if vartype.lower() == 'float':
            return float(vars(RawClass)[varname][0])
        if vartype.lower() == 'int_list':
            return np.asarray([int(i) for i in vars(RawClass)[varname]])
        if vartype.lower() == 'float_list':
            return np.asarray([float(i) for i in vars(RawClass)[varname]])
    except KeyError or AttributeError:
        return '! Variable '+varname+' not found in input file'


#Routine to extract and sort into a structure the most common VMEC variables.
#===============================================================================
class ReadInputVMEC:
    def __init__(self, File):

        self.Raw = RawData(File)
        
        self.Control  = Control(self.Raw)
        self.Grid     = Grid(self.Raw)
        self.FreeB    = FreeB(self.Raw)
        self.Pressure = Pressure(self.Raw)
        self.Flow     = Flow(self.Raw)
        self.Current  = Current(self.Raw)
        self.Boundary = Boundary(self.Raw)
        

    def WriteInput(self,name):

        f = open(name, "w")

        f.write(' &INDATA \n')

        #CONTROL PARAMETERS
        f.write(' \n')
        f.write(' !CONTROL PARAMETERS \n !---------------------------- \n')
        WriteVar(f, 'PRECON_TYPE', self.Control.PRECON_TYPE, 'str')
        WriteVar(f, 'PREC2D_THRESHOLD', self.Control.PREC2D_THRESHOLD, 'float')
        WriteVar(f, 'DELT', self.Control.DELT, 'float')
        WriteVar(f, 'NSTEP', self.Control.NSTEP, 'int')
        WriteVar(f, 'NITER_ARRAY', self.Control.NITER_ARRAY, 'int_list')
        WriteVar(f, 'NS_ARRAY', self.Control.NS_ARRAY, 'int_list')
        WriteVar(f, 'FTOL_ARRAY', self.Control.FTOL_ARRAY, 'float_list')

        #GRID PARAMETERS
        f.write(' \n')
        f.write(' !GRID PARAMETERS \n !---------------------------- \n')
        WriteVar(f, 'LASYM', self.Grid.LASYM, 'str')
        WriteVar(f, 'LRFP', self.Grid.LRFP, 'str')
        WriteVar(f, 'NFP', self.Grid.NFP, 'int')
        WriteVar(f, 'MPOL', self.Grid.MPOL, 'int')
        WriteVar(f, 'NTOR', self.Grid.NTOR, 'int')
        WriteVar(f, 'NTHETA', self.Grid.NTHETA, 'int')
        WriteVar(f, 'NZETA', self.Grid.NZETA, 'int')

        #FREE BOUNDARY PARAMETERS
        f.write(' \n')
        f.write(' !FREE BOUNDARY PARAMETERS \n !---------------------------- \n')
        WriteVar(f, 'LFREEB', self.FreeB.LFREEB, 'str')
        WriteVar(f, 'MGRID_FILE', self.FreeB.MGRID_FILE, 'str')
    
        WriteExcur(f, self.FreeB.EXTCURc, self.FreeB.EXTCUR)
            
        #PRESSURE PARAMETERS
        f.write(' \n')
        f.write(' !PRESSURE PARAMETERS \n !---------------------------- \n')
        WriteVar(f, 'GAMMA', self.Pressure.GAMMA, 'long_float')
        WriteVar(f, 'PRES_SCALE', self.Pressure.PRES_SCALE, 'str')
        WriteVar(f, 'PMASS_TYPE', self.Pressure.PMASS_TYPE, 'str')
        WriteVar(f, 'AM', self.Pressure.AM, 'long_float_list')
        WriteVar(f, 'AM_AUX_S', self.Pressure.AM_AUX_S, 'long_float_list')
        WriteVar(f, 'AM_AUX_F', self.Pressure.AM_AUX_F, 'long_float_list')
        
        #FLOWPARAMETERS
        f.write(' \n')
        f.write(' !FLOW PARAMETERS \n !---------------------------- \n')
        WriteVar(f, 'AT', self.Flow.AT, 'long_float_list')
        WriteVar(f, 'AH', self.Flow.AH, 'long_float_list')
        WriteVar(f, 'BCRIT', self.Flow.bcrit, 'long_float')
        
        
        #CURRENT PARAMETERS
        f.write(' \n')
        f.write(' !CURRENT PARAMETERS \n !---------------------------- \n')
        WriteVar(f, 'NCURR', self.Current.NCURR, 'int')
        WriteVar(f, 'AI', self.Current.AI, 'long_float_list')
        WriteVar(f, 'AC', self.Current.AC, 'long_float_list')
        WriteVar(f, 'CURTOR', self.Current.CURTOR, 'long_float')
        WriteVar(f, 'PIOTA_TYPE', self.Current.PIOTA_TYPE, 'str')
        WriteVar(f, 'PCURR_TYPE', self.Current.PCURR_TYPE, 'str')
        WriteVar(f, 'AI_AUX_S', self.Current.AI_AUX_S, 'long_float_list')
        WriteVar(f, 'AI_AUX_F', self.Current.AI_AUX_F, 'long_float_list')
        WriteVar(f, 'AC_AUX_S', self.Current.AC_AUX_S, 'long_float_list')
        WriteVar(f, 'AC_AUX_F', self.Current.AC_AUX_F, 'long_float_list')
        
        #BOUNDARY PARAMETERS
        f.write(' \n')
        f.write(' !BOUNDARY PARAMETERS \n !---------------------------- \n')
        WriteVar(f, 'PHIEDGE', self.Boundary.PHIEDGE, 'long_float')
        WriteVar(f, 'RAXIS', self.Boundary.RAXIS, 'long_float_list')
        WriteVar(f, 'ZAXIS', self.Boundary.ZAXIS, 'long_float_list')
        
        WriteBound(f, self.Boundary.RBCn, self.Boundary.RBCm, self.Boundary.RBC, 'RBC')
        WriteBound(f, self.Boundary.ZBSn, self.Boundary.ZBSm, self.Boundary.ZBS, 'ZBS')
        WriteBound(f, self.Boundary.RBSn, self.Boundary.RBSm, self.Boundary.RBS, 'RBS')
        WriteBound(f, self.Boundary.ZBCn, self.Boundary.ZBCm, self.Boundary.ZBC, 'ZBC')
        
        #End line
        f.write( '\n / \n &END \n')

        f.close()








        
class Control:
    def __init__(self,Raw):
        self.PRECON_TYPE      = GetVar(Raw,'PRECON_TYPE','str')
        self.PREC2D_THRESHOLD = GetVar(Raw,'PREC2D_THRESHOLD','float')
        self.DELT             = GetVar(Raw,'DELT','float')
        self.NSTEP            = GetVar(Raw,'NSTEP','int')
        self.NITER_ARRAY      = GetVar(Raw,'NITER_ARRAY','int_list')
        self.NS_ARRAY         = GetVar(Raw,'NS_ARRAY','int_list')
        self.FTOL_ARRAY       = GetVar(Raw,'FTOL_ARRAY','float_list')

class Grid:
    def __init__(self,Raw):
        self.LASYM  = GetVar(Raw,'LASYM','str')
        self.LRFP   = GetVar(Raw,'LRFP','str')
        self.NFP    = GetVar(Raw,'NFP','int')
        self.MPOL   = GetVar(Raw,'MPOL','int')
        self.NTOR   = GetVar(Raw,'NTOR','int')
        self.NTHETA = GetVar(Raw,'NTHETA','int')
        self.NZETA  = GetVar(Raw,'NZETA','int')

class FreeB:
    def __init__(self,Raw):
        self.LFREEB     = GetVar(Raw,'LFREEB','str')
        self.MGRID_FILE = GetVar(Raw,'MGRID_FILE','str')
        
        #Exception
        self.EXTCUR  = GetVar(Raw,'EXTCUR','float_list')
        self.EXTCURc = GetVar(Raw,'EXTCURc','int_list')

class Pressure:
    def __init__(self,Raw):
        self.GAMMA      = GetVar(Raw,'GAMMA','float')
        self.PRES_SCALE = GetVar(Raw, 'PRES_SCALE','float')
        self.PMASS_TYPE = GetVar(Raw,'PMASS_TYPE','str')
        self.AM         = GetVar(Raw,'AM','float_list')
        self.AM_AUX_S   = GetVar(Raw,'AM_AUX_S','float_list')
        self.AM_AUX_F   = GetVar(Raw,'AM_AUX_F','float_list')

class Flow:
	def __init__(self,Raw):
		self.AT			= GetVar(Raw,'AT','float_list')
		self.AH			= GetVar(Raw,'AH','float_list')
		self.bcrit		= GetVar(Raw,'BCRIT','float')
        
class Current:
    def __init__(self,Raw):
        self.NCURR      = GetVar(Raw,'NCURR','int')
        self.AI         = GetVar(Raw,'AI','float_list')
        self.AC         = GetVar(Raw,'AC','float_list')
        self.CURTOR     = GetVar(Raw,'CURTOR','float')
        self.PIOTA_TYPE = GetVar(Raw,'PIOTA_TYPE','str')
        self.PCURR_TYPE = GetVar(Raw,'PCURR_TYPE','str')
        self.AI_AUX_S   = GetVar(Raw,'AI_AUX_S','float_list')
        self.AI_AUX_F   = GetVar(Raw,'AI_AUX_F','float_list')
        self.AC_AUX_S   = GetVar(Raw,'AC_AUX_S','float_list')
        self.AC_AUX_F   = GetVar(Raw,'AC_AUX_F','float_list')
        
class Boundary:
    def __init__(self,Raw):
        self.PHIEDGE = GetVar(Raw,'PHIEDGE','float')
        self.RAXIS = GetVar(Raw,'RAXIS','float_list')
        self.ZAXIS = GetVar(Raw,'ZAXIS','float_list')
        
        #Exception
        self.RBCn, self.RBCm, self.RBC = GetVar(Raw,'RBCn','float_list'), GetVar(Raw,'RBCm','float_list'), GetVar(Raw,'RBC','float_list')
        self.RBSn, self.RBSm, self.RBS = GetVar(Raw,'RBSn','float_list'), GetVar(Raw,'RBSm','float_list'), GetVar(Raw,'RBS','float_list')
        self.ZBCn, self.ZBCm, self.ZBC = GetVar(Raw,'ZBCn','float_list'), GetVar(Raw,'ZBCm','float_list'), GetVar(Raw,'ZBC','float_list')
        self.ZBSn, self.ZBSm, self.ZBS = GetVar(Raw,'ZBSn','float_list'), GetVar(Raw,'ZBSm','float_list'), GetVar(Raw,'ZBS','float_list')
        

#*****************************************************************************
#**                        WRITE VMEC INPUT FILE                            **  
#**                                                                         **
#** This script writes an input file for VMEC provided all input parameters **
#** are specified.                                                          **
#**                                                                         **
#** AUTHOR: Guillermo Bustos Ramirez                                        ** 
#** LAST CHANGE: 23/09/2019                                                 **
#*****************************************************************************



def WriteVar(File,name, var, vartype):

    try:
        dummy = var[0]
    except TypeError:
        dummy = 'dummystr'
    except IndexError:
        dummy = 'dummystr'
        
    if dummy == '!':
        return None
    else:
        if vartype.lower() == 'str':
            File.write(' %s = %s , \n' %(name,var))
        if vartype.lower() == 'int':
            File.write(' %s = %i , \n' %(name,var))
        if vartype.lower() == 'float':
            File.write(' %s = %.2E , \n' %(name,var))
        if vartype.lower() == 'long_float':
            File.write(' %s = %.16E , \n' %(name,var))
        if vartype.lower() == 'int_list':
            File.write(' '+name+' = ')
            np.savetxt(File, var, fmt='%i', delimiter=' ', newline=" ")
            File.write(', \n')
        if vartype.lower() == 'float_list':
            File.write(' '+name+' = ')
            np.savetxt(File, var, fmt='%.2E', delimiter=' ', newline=" ")
            File.write(', \n')
        if vartype.lower() == 'long_float_list':
            File.write(' '+name+' = ')
            np.savetxt(File, var, fmt='%.16E', delimiter=' ', newline=" ")
            File.write(', \n')
            
def WriteBound(File, N, M, var, name):
    
    try:
        dummy = var[0]
    except TypeError:
        dummy = 'dummystr'
    
    if dummy == '!':
        return None
    else:
        for i in range(len(var)):
            File.write(' %s(%i,%i) = %.16E , \n' %(name, N[i], M[i], var[i]))
            
def WriteExcur(File, Ncur, EXTCUR):
    try:
        dummy = EXTCUR[0]
    except TypeError:
        dummy = 'dummystr'
    
    if dummy == '!':
        return None
    else:
        for i,Current in enumerate(EXTCUR):
            File.write(' EXTCUR(%i) = %.4E ,\n' %(Ncur[i],Current))





            
#===============================================================================






    


import numpy as np
import time
import sys
sys.path.insert(0, '../bin/')
import library as lib

'''
    Some tests to check the accuracy of Gauss-Legendre and Gauss-Jacobi methods of ~1/s functions near and at the axis.
'''


class GaussJacobi:
    def __init__(self,p):
        
        #Select the positiosn of the Gauss nodes.
        if p == 2: self.Pcoeff = [-0.999995000006250,0.333337222198843]
        if p == 3: self.Pcoeff = [-0.999997777775309,-0.289894139677351,0.689899117438264]
        if p == 4: self.Pcoeff = [-0.999998749997266,-0.575316121535975,0.181068260268316,0.822824570439263]
        if p == 5: self.Pcoeff = [-0.999999199997760,-0.720478227479051,-0.167178894832870,0.446315083462093,0.885791856124685]
    
        self.Weight = np.ones(shape=(p,1))	
        if p == 2: self.Weight[:,0] = [99999.5682,1.12499420424522]
        if p == 3: self.Weight[:,0] = [99998.804,1.44339981684823,0.44547474733231]
        if p == 4: self.Weight[:,0] = [99998.245,1.54864169207012,0.65735975042087,0.24189156392928]
        if p == 5: self.Weight[:,0] = [99997.807,1.5963072438726,0.74884187938281,0.38906675121166,0.15241784520583]

def f(x):
	return 0.5/x

def g(x):
	return 0.5 

NGauss = 5
GaussL = lib.GaussQuadrature(NGauss)
GaussJ = GaussJacobi(NGauss)



GJ = 0.
alpha = 0.
beta = -1.+1.E-05

GL = 0.
x1 = 0
x2 = 0.0196078431372549
h = x2-x1
for i in range(NGauss):
	node = x1 + 0.5*h*(1.+GaussL.Pcoeff[i])
	GL += 0.5*h*GaussL.Weight[i]*f(node)

	node = x1 + 0.5*h*(1.+GaussJ.Pcoeff[i])
	GJ += (0.5*h)**(alpha+beta+1.)*GaussJ.Weight[i]*g(node)

print ('GaussLegendre = %.10f'%(GL))
print ('GaussJacobi   = %.10f'%(GJ))







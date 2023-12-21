from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../bin/')
import time
import library as lib
from scipy import interpolate

a = 2000.
x  = np.linspace(1./(5.*a),1.,int(a))
y = 1./x

interp1 = interpolate.interp1d(x,y,kind='linear')
interp2 = interpolate.interp1d(x,y,kind='quadratic')
interp3 = interpolate.interp1d(x,y,kind='cubic')

x_ = np.linspace(1./(5.*a),1.,int(5.*a))
y1 = interp1(x_)
y2 = interp2(x_)
y3 = interp3(x_)

plt.plot(x,y, 'x-', label='Data')
#plt.plot(x_,y1, 'o-', label='lienar')
#plt.plot(x_,y2, '>-', label='quadratic')
plt.plot(x_,y3, '>-', label='cubic')
plt.grid()
plt.legend()

plt.show()

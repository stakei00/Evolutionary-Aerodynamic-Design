from xfoil import XFoil
from xfoil import model
import numpy as np 
import airfoil as af

xf = XFoil()
from xfoil.test import naca0012
#xf.airfoil = naca0012
#x = np.array([1.0000,0.9500,0.9000,0.8000,0.7000,0.6000,0.5000,0.4000,0.3000,0.2500,0.2000,0.1500,0.1000,0.0750,0.0500,0.0250,0.0125,0.0000,0.0125,0.0250,0.0500,0.0750,0.1000,0.1500,0.2000,0.2500,0.3000,0.4000,0.5000,0.6000,0.7000,0.8000,0.9000,0.9500,1.0000])
#y = np.array([0.0013,0.0114,0.0208,0.0375,0.0518,0.0636,0.0724,0.0780,0.0788,0.0767,0.0726,0.0661,0.0563,0.0496,0.0413,0.0299,0.0215,0.0000,-0.016,-0.022,-0.030,-0.034,-0.037,-0.041,-0.042,-0.042,-0.041,-0.038,-0.033,-0.027,-0.021,-0.015,-0.008,-0.004,-0.0013])

n2412 = af.Airfoil(None)
n2412.get_surface_coords(m=0.02, p=0.4, t=0.12, nPts=200)

naca2412 = model.Airfoil(x=n2412.x,y=n2412.y)

xf.airfoil = naca2412
xf.Re = 1e6
xf.max_iter = 40 
a, cl, cd, cm, co = xf.aseq(-20, 20, .5)

import matplotlib.pyplot as plt 
fig = plt.figure()
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

ax1.plot(a, cl)
ax1.set_xlabel("aoa"), ax1.set_ylabel("C_l")
ax2.plot(cl, cd)
ax2.set_xlabel("C_l"), ax2.set_ylabel("C_d")
ax3.plot(a, cm)
ax3.set_xlabel("aoa"), ax3.set_ylabel("C_m")

plt.show() 

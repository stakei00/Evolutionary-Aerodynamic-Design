import aero_evo as evo
from xfoil import XFoil

xf = XFoil()
from xfoil.test import naca0012
xf.airfoil = naca0012
xf.Re = 1e6
xf.max_iter = 40 
a, cl, cd, cm, co = xf.aseq(-20, 20, 1)

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

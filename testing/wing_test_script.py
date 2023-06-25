import sys, os
sys.path.append(os.getcwd())
import airfoil as af
import wing

"""
purpose of this script is to create a wing using two airfoils and analyze using 
XFOIL and AVL
"""

afoil_root = af.Airfoil(m=0.06, p=0.4, t=0.12, Re=6e5)
afoil_tip = af.Airfoil(m=0.02, p=0.2, t=0.14, Re=6e5)

import matplotlib.pyplot as plt 
plt.figure()
plt.plot(afoil_root.cd_interp, afoil_root.cl_interp, label="root")
plt.scatter(afoil_root.drag_coefficient, afoil_root.lift_coefficient)
plt.plot(afoil_tip.cd_interp, afoil_tip.cl_interp, label="tip")
plt.scatter(afoil_tip.drag_coefficient, afoil_tip.lift_coefficient)
plt.grid(), plt.legend()
plt.show()

taper = 0.5
AR = 12
sweep_deg = 0
twist_deg = -5

test_wing = wing.Wing(afoil_root, afoil_tip, taper, AR, sweep_deg, twist_deg)



test_wing.plot_avl_xfoil_results() 
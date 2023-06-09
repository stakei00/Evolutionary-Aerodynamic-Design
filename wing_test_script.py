import airfoil as af
import wing
import matplotlib.pyplot as plt 

"""
purpose of this scripy is to create a wing using two airfoils and analyze using 
XFOIL and AVL
"""

afoil_root = af.Airfoil(m=0.02, p=0.4, t=0.12, Re=3e5)
afoil_tip = af.Airfoil(m=0.04, p=0.4, t=0.12, Re=3e5)
span = 10
taper = 0.5
AR = 12
sweep_deg = 10
twist_deg = -5

test_wing = wing.Wing(afoil_root, afoil_tip, span, taper, AR, sweep_deg, twist_deg)
"""
fig,axs = plt.subplots(2,1, figsize=(16,6))
test_wing.plot_wing_planform(axs[0], linecolor="black")
test_wing.plot_wing_airfoils(axs[1], linecolor="black")

fig, axs = plt.subplots(2,1)
axs[0].plot(test_wing.alpha, test_wing.lift_coefficient)
axs[1].plot(test_wing.drag_coefficient, test_wing.lift_coefficient)
axs[0].set_xlabel("alpha"), axs[0].set_ylabel("CL")
axs[1].set_xlabel("CD"), axs[1].set_ylabel("CL")
plt.show()
"""
test_wing.plot_avl_xfoil_results()
pass 
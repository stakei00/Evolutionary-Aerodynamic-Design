import airfoil as af
import wing

"""
purpose of this scripy is to create a wing using two airfoils and analyze using 
XFOIL and AVL
"""

afoil_root = af.Airfoil(m=0.02, p=0.4, t=0.12)
afoil_tip = af.Airfoil(m=0.04, p=0.4, t=0.12)
span = 10
taper = 0.5
AR = 12
sweep_deg = 10
twist_deg = -5

test_wing = wing.Wing(afoil_root, afoil_tip, span, taper, AR, sweep_deg, twist_deg)

pass 
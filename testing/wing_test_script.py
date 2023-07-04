import sys, os
sys.path.append(os.getcwd())
import airfoil as af
import wing

"""
purpose of this script is to create a wing using two airfoils and analyze using 
XFOIL and AVL
"""

avl_settings = {
        "Nchord": 12, 
        "Nspan": 20, 
        "Sspace": 1, 
        "Cspace": 1, 
        "alpha_i": -10.0, 
        "alpha_f": 23.0, 
        "alpha_step": 2.0
    }

afoil_root = af.Airfoil(m=0.04, p=0.4, t=0.12, Re=5e5, search_airfoils=False)
afoil_tip = af.Airfoil(m=0.02, p=0.4, t=0.12, Re=5e5, search_airfoils=False)


taper = 0.5
AR = 12
sweep_deg = 0
twist_deg = 0

test_wing = wing.Wing(afoil_root, afoil_tip, taper, AR, sweep_deg, twist_deg, 
                      avl_settings=avl_settings)

test_wing.plot_avl_xfoil_results() 
import unittest
import sys, os
sys.path.append(os.getcwd())
import airfoil
import numpy as np
"""
tests airfoil class methods for a few well-known airfoils
"""

class Test_Airfoil(unittest.TestCase):

    def test_naca1412(self):
        print("testing naca1412")
        #test a few different reynolds numbers
        n1412_1e6 = airfoil.Airfoil(m=0.01, p=0.4, t=0.12, Re=1e6)
        n1412_5e5 = airfoil.Airfoil(m=0.01, p=0.4, t=0.12, Re=5e5)
        n1412_1e5 = airfoil.Airfoil(m=0.01, p=0.4, t=0.12, Re=1e5)

        #get max cl_cd of each 
        max_clcd_1e6 = max(np.divide(n1412_1e6.lift_coefficient, n1412_1e6.drag_coefficient))
        max_clcd_5e5 = max(np.divide(n1412_5e5.lift_coefficient, n1412_5e5.drag_coefficient))
        max_clcd_1e5 = max(np.divide(n1412_1e5.lift_coefficient, n1412_1e5.drag_coefficient))

        self.assertAlmostEqual(max_clcd_1e6, 83.8, places=-1)
        self.assertAlmostEqual(max_clcd_5e5, 73.5, places=-1)
        self.assertAlmostEqual(max_clcd_1e5, 44.1, places=-1)

    def test_naca2412(self):
        print("testing naca2412")
        #test a few different reynolds numbers
        n2412_1e6 = airfoil.Airfoil(m=0.02, p=0.4, t=0.12, Re=1e6)
        n2412_5e5 = airfoil.Airfoil(m=0.02, p=0.4, t=0.12, Re=5e5)
        n2412_1e5 = airfoil.Airfoil(m=0.02, p=0.4, t=0.12, Re=1e5)

        #get max cl_cd of each 
        max_clcd_1e6 = max(np.divide(n2412_1e6.lift_coefficient, n2412_1e6.drag_coefficient))
        max_clcd_5e5 = max(np.divide(n2412_5e5.lift_coefficient, n2412_5e5.drag_coefficient))
        max_clcd_1e5 = max(np.divide(n2412_1e5.lift_coefficient, n2412_1e5.drag_coefficient))

        self.assertAlmostEqual(max_clcd_1e6, 101.4, places=-1)
        self.assertAlmostEqual(max_clcd_5e5, 87.3, places=-1)
        self.assertAlmostEqual(max_clcd_1e5, 50, places=-1)
    
    def test_naca4412(self):
        print("testing naca4412")
        #test a few different reynolds numbers 
        n4412_1e6 = airfoil.Airfoil(m=0.04, p=0.4, t=0.12, Re=1e6)
        n4412_5e5 = airfoil.Airfoil(m=0.04, p=0.4, t=0.12, Re=5e5, aseq=[-12,15,1])
        n4412_1e5 = airfoil.Airfoil(m=0.04, p=0.4, t=0.12, Re=1e5)

        #get max cl_cd of each 
        max_clcd_1e6 = max(np.divide(n4412_1e6.lift_coefficient, n4412_1e6.drag_coefficient))
        max_clcd_5e5 = max(np.divide(n4412_5e5.lift_coefficient, n4412_5e5.drag_coefficient))
        max_clcd_1e5 = max(np.divide(n4412_1e5.lift_coefficient, n4412_1e5.drag_coefficient))

        self.assertAlmostEqual(max_clcd_1e6, 129.4, places=-1)
        self.assertAlmostEqual(max_clcd_5e5, 107.5, places=-1)
        self.assertAlmostEqual(max_clcd_1e5, 56.1, places=-1)
    
    def test_naca6412(self):
        print("testing naca6412")
        #test a few different reynolds numbers 
        n6412_1e6 = airfoil.Airfoil(m=0.06, p=0.4, t=0.12, Re=1e6, aseq=[-10,15,1])
        n6412_5e5 = airfoil.Airfoil(m=0.06, p=0.4, t=0.12, Re=5e5, aseq=[-10,15,1])
        n6412_1e5 = airfoil.Airfoil(m=0.06, p=0.4, t=0.12, Re=1e5, aseq=[-10,15,1])

        #get max cl_cd of each 
        max_clcd_1e6 = max(np.divide(n6412_1e6.lift_coefficient, n6412_1e6.drag_coefficient))
        max_clcd_5e5 = max(np.divide(n6412_5e5.lift_coefficient, n6412_5e5.drag_coefficient))
        max_clcd_1e5 = max(np.divide(n6412_1e5.lift_coefficient, n6412_1e5.drag_coefficient))

        self.assertAlmostEqual(max_clcd_1e6, 142.7, places=-1)
        self.assertAlmostEqual(max_clcd_5e5, 114.2, places=-1)
        self.assertAlmostEqual(max_clcd_1e5, 53.1, places=-1)
    
if __name__ == "__main__": 
    unittest.main()
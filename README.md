# Description
This is a simple optimization tool for designing a 2-section trapezoidal wing. It 
uses a genetic algorithm based on a scalar fitness/objective function. The main 
loop runs instances of XFOIL for the airfoil viscous-based drag and AVL for the 
lift-induced drag. The two airfoils are simple 3-parameter sections (NACA 4 digit).

# Requirements
Python 3.9
Packages: Numpy, Matplotlib, ... 
Custom Module: avlinput from OpenVSP with modified avlinput.py file 

# How to use
Intended use: create a python script, import aero_evo, build parameter 
dictionaries and call the function aero_evo.optimize(). Run time benefits from 
multiple cores. 
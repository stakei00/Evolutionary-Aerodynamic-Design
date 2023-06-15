import os, sys
sys.path.append(os.getcwd())
import aero_evo as evo
import json 

"""
list of commands for GUI widgets to be linked to
"""

def run_optimizer(): 
    """
    runs evolutionary optimizer on wing
    """
    settings = json.load(open("tkinter/user_settings.json", "r"))

    wing_parameters = settings["wing parameters"]
    study_parameters = settings["optimizer settings"]
    seed_wing = None
    
    def fitness(wing:object):
        """
        returns a fitness based on the highest lift-to-drag ratio
        """
        import numpy as np 
        lift_to_drag = np.divide(wing.lift_coefficient,wing.drag_coefficient)
        return max(lift_to_drag)
    
    evo.optimize(wing_parameters, study_parameters, fitness, seed_wing, live_plot=False)


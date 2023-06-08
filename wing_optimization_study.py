import aero_evo as evo 
"""

"""

wing_parameters = {
    "span":                     5, 
    "aspect ratio":             12,
    "taper":                    [0.2, 1],
    "sweep deg":                0,
    "twist deg":                [-10, 5],
    "root camber":              [0, 0.08],
    "root camber location":     [0.2, 0.6],
    "root thickness":           [0.04, 0.3],
    "tip camber":               [0, 0.08],
    "tip camber location":      [0.2, 0.6],
    "tip thickness":            [0.04, 0.3]
}

study_parameters = {
    "population size":              4,
    "children per generation":      2, 
    "gene mutation probability":    0.2,
    "child mutation probability":   0.4,
    "number of gens":               100
}

def fitness(wing:object):
    """
    returns a fitness based on the highest lift-to-drag ratio
    """
    import numpy as np 
    lift_to_drag = np.divide(wing.lift_coefficient,wing.drag_coefficient)
    return max(lift_to_drag)

optimized_wing = evo.optimize(wing_parameters=wing_parameters, \
                study_parameters=study_parameters, fitness_function=fitness, \
                    airfoil_history_path="airfoil_history.json", live_plot=True)
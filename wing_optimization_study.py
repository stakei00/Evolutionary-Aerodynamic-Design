import aero_evo as evo 
"""
script to run an optimization. Specifiy wing and study parameters, provide 
a seed wing (optional) and define your fitness function. The algorithm will 
MAXIMIZE the scalar fitness score, so if you are trying to minimize the scalar,  
it should return its inverse or negative. Fitness function could also 
return a generic score based on multiple desirable traits each with scoring 
weights. (the sky's the limit)
"""

#get reynolds number at design flight conditions (based on L_ref = 1)
altitude_ft = 1000
true_airspeed_ft_s = 80
Re = evo.get_reynolds_number(altitude_ft, true_airspeed_ft_s, units="US")

wing_parameters = {
    #defines the fixed and varying parameters. varying parameters must be 
    # specified as list 
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
    "tip thickness":            [0.04, 0.3],
    "reynolds number":          Re
}

study_parameters = {
    #defines genetic algorithm study settings
    "population size":              6,
    "children per generation":      3, 
    "gene mutation probability":    0.125,
    "child mutation probability":   1,
    "number of gens":               5
}

seed_wing = {
    #specifies traits of a seed wing to intialize population with. Optional input 
    "span":                     5, 
    "aspect ratio":             12,
    "taper":                    0.6,
    "sweep deg":                0,
    "twist deg":                -3,
    "root camber":              0.02,
    "root camber location":     0.4,
    "root thickness":           0.12,
    "tip camber":               0.02,
    "tip camber location":      0.4,
    "tip thickness":            0.12,
    "reynolds number":          Re
}

#user defined fitness function: 
def fitness(wing:object):
    """
    returns a fitness based on the highest lift-to-drag ratio
    """
    import numpy as np 
    lift_to_drag = np.divide(wing.lift_coefficient,wing.drag_coefficient)
    return max(lift_to_drag)

#optimization function call 
optimized_wing = evo.optimize(wing_parameters=wing_parameters, \
                    study_parameters=study_parameters,\
                        fitness_function=fitness,\
                            seed_wing=seed_wing,\
                                airfoil_history_path="airfoil_history.json",\
                                    live_plot=True)
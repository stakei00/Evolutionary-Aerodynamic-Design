import sys, os 
sys.path.append(os.getcwd())
import aero_evo as evo 
"""
script to run an optimization. Specifiy wing and study parameters, provide 
a seed wing (optional) and define your fitness function. The algorithm will 
MAXIMIZE the scalar fitness score, so if you are trying to minimize the scalar,  
fitness func should return its inverse or negative. Fitness function could also 
return a generic score based on multiple desirable traits each with scoring 
weights. (the sky's the limit)
"""

#get reynolds number to chord ratio at design flight conditions
altitude_ft = 1000
true_airspeed_ft_s = 80
Re = evo.get_reynolds_number(altitude_ft, true_airspeed_ft_s, units="US")

wing_parameters = {
    #defines the fixed and varying parameters. varying parameters must be 
    # specified as 2-element lists 
    "span":                     5,              #full wing span 
    "aspect ratio":             12,             #full aspect ratio
    "taper":                    [0.3, 1],       #taper ratio 
    "sweep deg":                0,             #wing sweep at qua
    "twist deg":                [-5, 5],       #wing tip incidenc
    "root camber":              [0, 0.06],      #NACA X...
    "root camber location":     [0.2, 0.6],     #NACA .X..
    "root thickness":           [0.08, 0.2],    #NACA ..XX
    "tip camber":               [0, 0.06],
    "tip camber location":      [0.2, 0.6],
    "tip thickness":            [0.08, 0.3],
    "reynolds number":          Re              #reynolds # to chord ratio (be aware of this)
}  

study_parameters = {
    #defines genetic algorithm study settings
    "population size":              10, 
    "children per generation":      5, #number of new chromosomes per generation (< population size)
    "gene mutation probability":    0.125, #probability that a gene will mutate 
    "child mutation probability":   0.5, #probabilty that child will have mutation(s)
    "number of gens":               400 #number of generations/iterations for study 
}

seed_wing = {
    #specifies traits of a seed wing to intialize population with. Optional input 
    "span":                     5,      #full wing span 
    "aspect ratio":             12,     #full aspect ratio 
    "taper":                    0.6,    #taper ratio 
    "sweep deg":                0,     #wing sweep at quarter-chord line 
    "twist deg":                -3,     #wing tip incidence angle (deg)
    "root camber":              0.02,   #NACA X...
    "root camber location":     0.4,    #NACA .X..
    "root thickness":           0.12,   #NACA ..XX
    "tip camber":               0.02, 
    "tip camber location":      0.4,
    "tip thickness":            0.12, 
    "reynolds number":          Re #reynolds # to chord ratio (be aware of this)
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
                            seed_wing=None,\
                                live_plot=True)
optimized_wing.plot_avl_xfoil_results()
optimized_wing.export_wing()
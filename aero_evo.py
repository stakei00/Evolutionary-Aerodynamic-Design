import wing
import random
import airfoil
import numpy as np
import multiprocessing
import time

class Population:
    """
    creates and manipulates population
    """
    def __init__(self, size:int, wing_parameters:dict, seed_wing:dict=None,\
                 airfoil_hist:dict=None) -> None: 
        """
        initializes population of chromosomes
        Inputs: 
            size: size of the population 
            wing_parameters: values and ranges for wing independent parameters 
            seed_wing: if provided, population initializes with clones of this wing 
            airfoil_hist: if provided, store XFOIL solutions for future reference 
                (speeds things up)
        """
        self.chroms = []
        self.airfoil_hist = airfoil_hist
        self.best_chrom_fitness = []
        
        self.variables = {} #variables and limits
        self.constants = {}
        for key in wing_parameters.keys():
            if not isinstance(wing_parameters[key], list): 
                self.constants[key] = wing_parameters[key]
                continue

            self.variables[key] = wing_parameters[key]

        def unpack_parameters_dict(parameters_dict):
            root, tip, span, taper, AR, sweep, twist = [None,None,None],\
                    [None,None,None],None,None,None,None,None 

            chrom_header = {}
            for key in parameters_dict.keys(): 

                    if isinstance(parameters_dict[key], list or tuple): 
                        val = random.uniform(parameters_dict[key][0], parameters_dict[key][-1])
                        chrom_header[key] = val
                    else:
                        val = parameters_dict[key]

                    if key == "span":                   span = val
                    elif key == "aspect ratio":         AR = val 
                    elif key == "taper":                taper = val
                    elif key == "sweep deg":            sweep = val
                    elif key == "twist deg":            twist = val 
                    elif key == "root camber":          root[0] = val
                    elif key == "root camber location": root[1] = val
                    elif key == "root thickness":       root[2] = val
                    elif key == "tip camber":           tip[0] = val
                    elif key == "tip camber location":  tip[1] = val
                    elif key == "tip thickness":        tip[2] = val

            if None in [root[0], root[1], root[2], tip[0], tip[1], tip[2],\
                        span, taper, AR, sweep, twist]: 
                raise ValueError("some input left unspecified")
            
            return chrom_header, root, tip, span, taper, AR, sweep, twist
        
        parameters_dict = wing_parameters
        if seed_wing is not None: #if a seed wing is provided use instead of rng 
            
            seed_chrom_header = {}
            for key in parameters_dict.keys():
                if isinstance(wing_parameters[key], list):
                    seed_chrom_header[key] = seed_wing[key] 

            parameters_dict = seed_wing

        while len(self.chroms) < size:
            
            chrom_header, root, tip, span, taper, AR, sweep, twist = unpack_parameters_dict(parameters_dict)
            if seed_wing is not None: chrom_header = seed_chrom_header

            airfoil_root = airfoil.Airfoil(root[0], root[1], root[2], 10e6, airfoil_hist=self.airfoil_hist)
            if not hasattr(airfoil_root, "CDCL_avl"): continue
            self.airfoil_hist = airfoil_root.try_add_to_hist(self.airfoil_hist) 

            airfoil_tip = airfoil.Airfoil(tip[0], tip[1], tip[2], 10e6, airfoil_hist=self.airfoil_hist)
            if not hasattr(airfoil_tip, "CDCL_avl"): continue
            self.airfoil_hist = airfoil_tip.try_add_to_hist(self.airfoil_hist)

            chrom = Chromosome(airfoil_root, airfoil_tip, span, taper, AR, sweep, twist)
            chrom.chrom_header = chrom_header
            self.chroms.append(chrom)

    def selection(self, k) -> None: 
        """
        reduces number of chromosomes using truncation-selection based on 
        individual fitness scores
        k: # of members to delete
        """
        chroms_sorted = sorted(self.chroms, key=lambda c: c.fitness)
        self.chroms = chroms_sorted[k:]

    def crossover_and_mutation(self, k, p_gene_mut, p_child_mut) -> None: 
        """
        creates new chromosomes and mutates using interpolated mutations
        k: number of chromosomes to generate
        p_gene_mut: probability that a gene in a mutated chromosome will experience 
                    mutation 
        p_child_mut: probability that a given child will experience mutation 
        """

        def extract_inputs(chrom_header:dict, constants:dict):

            root, tip, span, taper, AR, sweep, twist = [None,None,None],\
                [None,None,None],None,None,None,None,None
            
            for key in list(chrom_header.keys()) + list(constants.keys()):
                if key in chrom_header.keys(): dict_=chrom_header 
                if key in constants.keys(): dict_=constants

                if key == "span":                   span = dict_[key]                
                if key == "aspect ratio":           AR = dict_[key]       
                if key == "taper":                  taper = dict_[key]              
                if key == "sweep deg":              sweep = dict_[key]           
                if key == "twist deg":              twist = dict_[key]           
                if key == "root camber":            root[0] = dict_[key]         
                if key == "root camber location":   root[1] = dict_[key]
                if key == "root thickness":         root[2] = dict_[key]      
                if key == "tip camber":             tip[0] = dict_[key]          
                if key == "tip camber location":    tip[1] = dict_[key] 
                if key == "tip thickness":          tip[2] = dict_[key]       

            return root, tip, span, taper, AR, sweep, twist

        new_chroms = []
        while len(new_chroms) < k: 

            i1 = random.choice(range(len(self.chroms)))
            i2 = i1
            while i2 == i1: 
                i2 = random.choice(range(len(self.chroms)))

            parent1, parent2 = self.chroms[i1], self.chroms[i2]
            child_chrom_header = {}
            for i,_ in enumerate(parent1.chrom_header):
                key = list(parent1.chrom_header.keys())[i]
                child_chrom_header[key] = 0.5*(parent1.chrom_header[key] + parent2.chrom_header[key]) #averaged 
            
            #apply mutations to child_chrom_header
            if random.uniform(0,1) <= p_child_mut:
                
                for key in child_chrom_header.keys(): 
                    
                    if random.uniform(0,1) >= p_gene_mut: continue 
                    mutation_factor = np.random.normal(loc=1)
                    child_chrom_header[key] *= mutation_factor
                    
                    upper_lim, lower_lim = self.variables[key][-1], self.variables[key][0]
                    if child_chrom_header[key] < lower_lim: 
                        child_chrom_header[key] = lower_lim
                    elif child_chrom_header[key] > upper_lim:
                        child_chrom_header[key] = upper_lim  

            root, tip, span, taper, AR, sweep, twist = extract_inputs(child_chrom_header, self.constants)

            airfoil_root = airfoil.Airfoil(root[0], root[1], root[2], 10e6, airfoil_hist=self.airfoil_hist)
            if not hasattr(airfoil_root, "CDCL_avl"): continue
            self.airfoil_hist = airfoil_root.try_add_to_hist(self.airfoil_hist)
             
            airfoil_tip = airfoil.Airfoil(tip[0], tip[1], tip[2], 10e6, airfoil_hist=self.airfoil_hist)
            if not hasattr(airfoil_tip, "CDCL_avl"): continue
            self.airfoil_hist = airfoil_tip.try_add_to_hist(self.airfoil_hist)
            
            new_chrom = Chromosome(airfoil_root, airfoil_tip, span, taper, AR, sweep, twist)
            new_chrom.chrom_header = child_chrom_header
            new_chroms.append(new_chrom)

        [self.chroms.append(chrom) for chrom in new_chroms] #add children to population 

    def get_best_chrom(self) -> None:
        """
        finds best fitness chromosome in the population
        """
        chroms_sorted = sorted(self.chroms, key=lambda c: c.fitness)
        self.best_chrom = chroms_sorted[-1]
        self.best_chrom_fitness.append(self.best_chrom.fitness) 


class Chromosome(wing.Wing):
    """
    creates a single chromosome 
    """
    def __init__(self, airfoil_root, airfoil_tip, span, taper_ratio, \
                         aspect_ratio, sweep_deg, twist_deg) -> None: 
        """
        creates individual wing
        """
        super().__init__(airfoil_root, airfoil_tip, span, taper_ratio, \
                         aspect_ratio, sweep_deg, twist_deg)

    def evaluate_fitness(self, fitness_func) -> None: 
        """
        scores individual chromosome based on provided fitness function
        """
        self.fitness = fitness_func(self)


def optimize(wing_parameters:dict, study_parameters:dict, fitness_function,\
             seed_wing:dict=None, airfoil_history_path:str=None,\
                live_plot=False) -> object:
    """
    Runs genetic algorithm study
    """
    t_start = time.time()/60
    if live_plot: 
        import matplotlib.pyplot as plt
        fig, axs = initialize_plot()
      
    #open airfoil history file
    airfoil_hist = None 
    if airfoil_history_path is not None: 
        import json
        airfoil_hist = json.load(open(airfoil_history_path, "r"))

    #intialize population of m individuals
    popSize = study_parameters["population size"]
    k = study_parameters["children per generation"]
    p_gene_mut, p_child_mut = study_parameters["gene mutation probability"],\
                                study_parameters["child mutation probability"]
    population = Population(size=popSize, wing_parameters=wing_parameters,\
                            seed_wing=seed_wing, airfoil_hist=airfoil_hist)  
    [chrom.evaluate_fitness(fitness_function) for chrom in population.chroms]
    population.get_best_chrom()
    nIter = study_parameters["number of gens"]
    
    prev_best = population.best_chrom
    n = 0  
    while n < nIter:

        if live_plot: update_plot(fig, axs, population, prev_best, "running...", time0=t_start)

        prev_best = population.best_chrom
        population.selection(k) #kill off worst performing chromosomes
        population.crossover_and_mutation(k, p_gene_mut, p_child_mut) #generate new children 
        [chrom.evaluate_fitness(fitness_function) for chrom in population.chroms] #evaluate fitness 
        population.get_best_chrom()
        n += 1  
    
    if live_plot: 
        update_plot(fig, axs, population, prev_best, "finished", time0=t_start)
        plt.ioff()
        plt.show()

    #save airfoil history  
    json.dump(population.airfoil_hist, open(airfoil_history_path, "w"))

    #return    
    return population.best_chrom


def initialize_plot():
        
        import matplotlib.pyplot as plt
        plt.ion()
        #plt.style.use("dark_background")
        fig = plt.figure(figsize=(16,8))
        gs = fig.add_gridspec(3,4)

        ax1 = fig.add_subplot(gs[0,:3])
        ax1.set_title("planform")

        ax2 = fig.add_subplot(gs[1,:3])
        ax2.set_title("section")

        ax3 = fig.add_subplot(gs[2,:])
        ax3.set_xlabel("generation")
        ax3.set_ylabel("best fitness")
        ax3.grid()

        ax4 = fig.add_subplot(gs[:2,3])
        ax4.set_xlim(0,1), ax4.set_ylim(0,1)
        

        return fig, [ax1, ax2, ax3, ax4]


def update_plot(fig, axs, population, prev_best_chrom, status_msg, time0=None,):
    ax1, ax2, ax3, ax4 = axs
    ax1.clear(), ax2.clear(), ax4.clear()

    best_wing = population.best_chrom
    prev_best_chrom.plot_wing_planform(ax1, linecolor="lightgrey")
    best_wing.plot_wing_planform(ax1, linecolor="black")
    prev_best_chrom.plot_wing_airfoils(ax2, linecolor="lightgrey")
    best_wing.plot_wing_airfoils(ax2, linecolor="black")
    
    n = len(population.best_chrom_fitness)
    ax3.plot(range(n),population.best_chrom_fitness, "-o", color="red", linewidth=0.5)
    
    text =  f"status: {status_msg}\n"

    if time0 is not None:
        text += f"time elapsed: {round(time.time()/60 - time0,2)} mins\n"

    text += f"generation number: {n}\n" +\
            f"best fitness: {round(population.best_chrom_fitness[-1],3)}\n\n" +\
            f"root airfoil: NACA{best_wing.airfoil_root.NACA_4series_desig}\n" +\
            f"tip airfoil:  NACA{best_wing.airfoil_tip.NACA_4series_desig}\n" +\
            f"aspect ratio: {round(best_wing.AR,2)}\n" +\
            f"span: {round(best_wing.span, 2)}\n" +\
            f"taper ratio: {round(best_wing.taper,3)}\n" +\
            f"q-chord sweep: {round(best_wing.sweep_deg,2)} deg\n" +\
            f"tip twist: {round(best_wing.tip_twist_deg,2)} deg\n"
    
    try: text += f"max lift coefficient: {round(best_wing.max_lift_coefficient, 2)}\n"
    except: pass 
                
    ax4.invert_yaxis()
    ax4.set_xticks([]), ax4.set_yticks([])
    ax4.text(0.05, 0.03, text, fontsize=10, ha="left", va="top")
    fig.canvas.draw()
    fig.canvas.flush_events()
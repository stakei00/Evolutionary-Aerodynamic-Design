import wing
import random
import airfoil
import numpy as np
import multiprocessing as mp
import math 
import time
import os 

class Population:
    """
    creates and manipulates population
    """
    def __init__(self, size:int, wing_parameters:dict, mutation_probs:list, 
                avl_settings:dict, xfoil_settings:dict, seed_wing:dict=None, 
                multiproc:int=1) -> None: 
        """
        initializes population of chromosomes
        Inputs: 
            size: size of the population 
            wing_parameters: values and ranges for wing independent parameters 
            seed_wing: if provided, population initializes with clones of this wing 

        """
        self.chroms = []
        self.best_chrom_fitness = []
        self.best_chrom_history = []
        self.p_gene_mut = mutation_probs[1]
        self.p_child_mut = mutation_probs[0]
        self.pop_size = size
        self.xfoil_settings = xfoil_settings
        self.avl_settings = avl_settings
        
        self.multiproc = multiproc
        if multiproc > mp.cpu_count(): self.multiproc = mp.cpu_count()
        
        #separate constants from variables 
        self.variables = {} #variables and limits
        self.constants = {}
        for key in wing_parameters.keys():
            if not isinstance(wing_parameters[key], list): 
                self.constants[key] = wing_parameters[key]
                continue
            self.variables[key] = wing_parameters[key]
        
        parameters_dict = wing_parameters
        if seed_wing is not None: #if a seed wing is provided use instead of rng 
            seed_chrom_header = {}
            for key in parameters_dict.keys():
                if isinstance(wing_parameters[key], list):
                    seed_chrom_header[key] = seed_wing[key] 

        while len(self.chroms) < self.pop_size:
            
            if self.multiproc > 1: 
                num_to_gen = size - len(self.chroms)
                if self.multiproc < num_to_gen:
                    num_to_gen = self.multiproc

                if seed_wing is not None: 
                    headers = [seed_chrom_header for _ in range(num_to_gen)]
                else: 
                    headers = [self.create_header(mutate=False) \
                               for _ in range(num_to_gen)]

                pool = mp.Pool(processes=num_to_gen)
                chroms = pool.map(self.generate_chromosome, headers)
                [self.chroms.append(c) for c in chroms if c is not None]

            else: 
                if seed_wing is not None: 
                    header = seed_chrom_header
                else: 
                    header = self.create_header(mutate=False)
                chrom = self.generate_chromosome(header)
                
                if chrom is not None: 
                    self.chroms.append(chrom)

    def extract_inputs(self, chrom_header:dict) -> list:
            
        root, tip, taper, AR, sweep, twist, re_c = [None,None,None],\
            [None,None,None],None,None,None,None,None
        
        for key in list(chrom_header.keys()) + list(self.constants.keys()):
            if key in chrom_header.keys(): dict_=chrom_header 
            if key in self.constants.keys(): dict_=self.constants
            if key == "aspect ratio":           AR = dict_[key]       
            if key == "taper":                  taper = dict_[key]              
            if key == "sweep deg":              sweep = dict_[key]           
            if key == "twist deg":              twist = dict_[key]
            if key == "Re/chord":               re_c = dict_[key]           
            if key == "root camber":            root[0] = dict_[key]         
            if key == "root camber location":   root[1] = dict_[key]
            if key == "root thickness":         root[2] = dict_[key]      
            if key == "tip camber":             tip[0] = dict_[key]          
            if key == "tip camber location":    tip[1] = dict_[key] 
            if key == "tip thickness":          tip[2] = dict_[key]       
        
        return root, tip, taper, AR, sweep, twist, re_c
    
    def create_header(self, parents:list=None, mutate:bool=True) -> dict:
        """
        creates a header dictionary which can then be used to create a chromosome
        """
        header = {}
        if parents is None: 
            #generate header with random traits     
            for key in self.variables.keys():
                var = self.variables[key] 
                header[key] = random.uniform(var[0], var[1])

        else: 
            #generate header by averaging parents traits
            parent1, parent2 = parents[0], parents[1]
            for i,_ in enumerate(parent1.chrom_header):
                key = list(parent1.chrom_header.keys())[i]
                header[key] = 0.5*(parent1.chrom_header[key] + parent2.chrom_header[key]) #averaged 

        if mutate: 
            if random.uniform(0,1) <= self.p_child_mut:
                for key in header.keys(): 
                
                    if random.uniform(0,1) >= self.p_gene_mut: continue 
                    mutation_factor = np.random.normal(loc=0, scale=1)
                    header[key] *= 1 + mutation_factor

                    upper_lim, lower_lim = self.variables[key][-1], self.variables[key][0]
                    if header[key] < lower_lim: 
                        header[key] = lower_lim
                    elif header[key] > upper_lim:
                        header[key] = upper_lim 

        return header           

    def generate_chromosome(self, child_chrom_header) -> object:
        """
        multiprocessing target function for creating new chromosomes
        """

        #extract from settings 
        aseq = [self.xfoil_settings["alpha_i"],self.xfoil_settings["alpha_f"],
                self.xfoil_settings["alpha_step"]]
        n_iter = self.xfoil_settings["iterations"]
        ncrit = self.xfoil_settings["Ncrit"]
        t_timeout = self.xfoil_settings["timeout limit"]
        
        try:
            root, tip, taper, AR, sweep, twist, re_c = \
                self.extract_inputs(child_chrom_header)
            root_chord = 1
            Re_root = round(re_c*root_chord, -5)
            airfoil_root = airfoil.Airfoil(root[0], root[1], root[2], Re_root, 
                                           aseq=aseq, n_iter=n_iter, ncrit=ncrit, 
                                           t_timeout=t_timeout)
            if not hasattr(airfoil_root, "alpha_raw"): return None

            tip_chord = taper*root_chord
            Re_tip = round(re_c*tip_chord, -5)
            airfoil_tip = airfoil.Airfoil(tip[0], tip[1], tip[2], Re_tip, 
                                          aseq=aseq, n_iter=n_iter, ncrit=ncrit, 
                                          t_timeout=t_timeout)
            if not hasattr(airfoil_tip, "alpha_raw"): return None

            new_chrom = Chromosome(airfoil_root, airfoil_tip, taper, AR, sweep, 
                                   twist, self.avl_settings)
            new_chrom.chrom_header = child_chrom_header
            return new_chrom
        
        except: 
            return None

    def selection(self, k) -> None: 
        """
        reduces number of chromosomes using truncation-selection based on 
        individual fitness scores
        k: # of members to delete
        """
        chroms_sorted = sorted(self.chroms, key=lambda c: c.fitness)
        self.chroms = chroms_sorted[k:]

    def crossover_and_mutation(self, k) -> None: 
        """
        creates new chromosomes and mutates using interpolated mutations
        k: number of chromosomes to generate
        p_gene_mut: probability that a gene in a mutated chromosome will experience 
                    mutation 
        p_child_mut: probability that a given child will experience mutation 
        """

        new_chroms = []
        while len(new_chroms) < k: 
            
            if self.multiproc > 1: 
                num_to_gen = k - len(new_chroms)
                if self.multiproc < num_to_gen:
                    num_to_gen = self.multiproc

                headers = []    
                for _ in range(num_to_gen):
                    i1 = random.choice(range(len(self.chroms)))
                    i2 = i1
                    while i2 == i1: 
                        i2 = random.choice(range(len(self.chroms)))
                    parents = [self.chroms[i1], self.chroms[i2]]
                    headers.append(self.create_header(parents))

                pool = mp.Pool(processes=num_to_gen)
                chroms = pool.map(self.generate_chromosome, headers)
                [new_chroms.append(c) for c in chroms if c is not None]

            else: 
                i1 = random.choice(range(len(self.chroms)))
                i2 = i1
                while i2 == i1: 
                    i2 = random.choice(range(len(self.chroms)))
                parents = [self.chroms[i1], self.chroms[i2]]
                header = self.create_header(parents)
                chrom = self.generate_chromosome(header)
                if chrom is not None: 
                    new_chroms.append(chrom)

        [self.chroms.append(chrom) for chrom in new_chroms] #add children to population 

    def get_best_chrom(self) -> None:
        """
        finds best fitness chromosome in the population
        """
        chroms_sorted = sorted(self.chroms, key=lambda c: c.fitness)
        self.best_chrom = chroms_sorted[-1]
        
        if self.best_chrom not in self.best_chrom_history:
            self.best_chrom_history.append(self.best_chrom)

        self.best_chrom_fitness.append(self.best_chrom.fitness) 


class Chromosome(wing.Wing):
    """
    creates a single chromosome 
    """
    def __init__(self, airfoil_root, airfoil_tip, taper_ratio, \
                         aspect_ratio, sweep_deg, twist_deg, avl_settings) -> None: 
        """
        creates individual wing
        """
        super().__init__(airfoil_root, airfoil_tip, taper_ratio, \
                         aspect_ratio, sweep_deg, twist_deg, avl_settings)

    def evaluate_fitness(self, fitness_func) -> None: 
        """
        scores individual chromosome based on provided fitness function
        """
        self.fitness = fitness_func(self)


def optimize(wing_parameters:dict, study_parameters:dict, fitness_function,\
             seed_wing:dict=None,live_plot=False, multiproc=True) -> object:
    """
    Runs genetic algorithm study
    """
    t_start = time.time()/60
    if live_plot: 
        import matplotlib.pyplot as plt
        fig, axs = initialize_3figplot()

    #intialize population of m individuals
    popSize = study_parameters["population size"]
    k = study_parameters["children per generation"]
    if k >= popSize: 
        raise ValueError("invalid children per generation and population size")
    
    p_gene_mut, p_child_mut = study_parameters["gene mutation probability"],\
                                study_parameters["child mutation probability"]
    
    #clear existing xfoil polar files
    file_list = os.listdir(os.getcwd())
    [os.remove(file) for file in file_list if file.split("_")[0] == "polar"]

    #create population 
    population = Population(size=popSize, wing_parameters=wing_parameters,\
                            mutation_probs=[p_child_mut, p_gene_mut],\
                                seed_wing=seed_wing, multiproc=6)  
    
    #evaluate fitness of initial population 
    [chrom.evaluate_fitness(fitness_function) for chrom in population.chroms]
    population.get_best_chrom()

    nIter = study_parameters["number of gens"]
    prev_best = population.best_chrom
    
    n = 0  
    while n < nIter:
        #run genetic algorithm:
        if live_plot: update_3figplot(fig, axs, population, prev_best, "running...",\
                                   time0=t_start)

        prev_best = population.best_chrom
        population.selection(k) #kill off worst performing chromosomes
        population.crossover_and_mutation(k) #generate new children 
        [chrom.evaluate_fitness(fitness_function) for chrom in population.chroms] #evaluate fitness 
        population.get_best_chrom()
        n += 1  

    if live_plot: 
        update_3figplot(fig, axs, population, prev_best, "finished", time0=t_start)
        plt.ioff()
        plt.show()

    #return    
    return population.best_chrom


def initialize_3figplot():
        
        import matplotlib
        import matplotlib.pyplot as plt
        #plt.ion()
        #plt.style.use("dark_background")
        matplotlib.rcParams["toolbar"] = "None" #hide the toolbar
        fig = plt.figure(figsize=(12,6), dpi=100)
        fig.canvas.manager.set_window_title("live plot")
        gs = fig.add_gridspec(3,4)

        ax1 = fig.add_subplot(gs[0,:3])
        ax1.spines[:].set_color("lightgrey")
         
        ax2 = fig.add_subplot(gs[1,:3])
        ax2.spines[:].set_color("lightgrey")

        ax3 = fig.add_subplot(gs[2,:])
        ax3.spines[:].set_color("lightgrey")
        ax3.set_xlabel("generations")
        ax3.set_ylabel("fitness")
        ax3.grid(color="lightgrey")

        ax4 = fig.add_subplot(gs[:2,3])
        ax4.set_xlim(0,1), ax4.set_ylim(0,1)
        ax4.spines[:].set_color("lightgrey")
                
        return fig, [ax1, ax2, ax3, ax4]


def update_3figplot(fig, axs, population, prev_best_chrom, status_msg, time0=None,):
    ax1, ax2, ax3, ax4 = axs
    ax1.clear(), ax2.clear(), ax4.clear()

    best_wing = population.best_chrom
    prev_best_chrom.plot_wing_planform(ax1, linecolor="dimgrey")
    best_wing.plot_wing_planform(ax1, linecolor="black")
    prev_best_chrom.plot_wing_airfoils(ax2, linecolor="dimgrey")
    best_wing.plot_wing_airfoils(ax2, linecolor="black")
    
    n = len(population.best_chrom_fitness)
    ax3.plot(range(n),population.best_chrom_fitness, "-o", color="red",\
             linewidth=0.5, markersize=3)
    
    text =  f"status: {status_msg}\n"

    if time0 is not None:
        text += f"time elapsed: {round(time.time()/60 - time0,2)} mins\n"

    text += f"generation number: {n}\n" +\
            f"best fitness: {round(population.best_chrom_fitness[-1],3)}\n\n" +\
            f"root airfoil: NACA{best_wing.airfoil_root.NACA_4series_desig}\n" +\
            f"tip airfoil:  NACA{best_wing.airfoil_tip.NACA_4series_desig}\n" +\
            f"aspect ratio: {round(best_wing.AR,2)}\n" +\
            f"span/root_chord: {round(best_wing.span, 2)}\n" +\
            f"Reynolds/chord: {round(best_wing.Re_c, 5)}\n" +\
            f"taper ratio: {round(best_wing.taper,3)}\n" +\
            f"q-chord sweep: {round(best_wing.sweep_deg,2)} deg\n" +\
            f"tip twist: {round(best_wing.tip_twist_deg,2)} deg\n"
    
    try: text += f"max lift coefficient: {round(best_wing.max_lift_coefficient, 2)}\n"
    except: pass 
                
    ax4.invert_yaxis()
    ax4.set_xticks([]), ax4.set_yticks([])
    ax4.text(0.05, 0.03, text, fontsize=10, ha="left", va="top", color="black")
    fig.canvas.draw()
    fig.canvas.flush_events()


def get_reynolds_number(h:float , V:float, units:str, l:float=1) -> float:
    
    def get_reynolds_number_SI(h_m, V_m_s, l_m):
        dens_kg_m3 = (1.04884 - 23.659414e-6*h_m)**4.2558797
        temp_k = 288.15 - 6.5e-3*h_m
        temp_ref_k = 273.15
        S_k = 110.4
        dynVsc_pas = 1.716e-5*(temp_k/temp_ref_k)**(3/2)*(temp_ref_k + S_k)/(temp_k + S_k)
        Re = dens_kg_m3*V_m_s*l_m/dynVsc_pas
        return Re   

    def get_reynolds_number_US(h_ft, V_ft_s, l_ft): 
        h_m = 0.3048*h_ft
        V_m_s = 0.3048*V_ft_s
        l_m = 0.3048*l_ft
        return get_reynolds_number_SI(h_m, V_m_s, l_m)

    if units in ["SI", "metric"]:
        return get_reynolds_number_SI(h, V, l)
    
    elif units in ["US", "imperial"]:
        return get_reynolds_number_US(h, V, l)
    else: 
        raise ValueError("unknown units")
    
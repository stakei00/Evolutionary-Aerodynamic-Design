from tkinter import *
from tkinter import scrolledtext
from tkinter.ttk import *
import json
import threading
import os, sys
import time
sys.path.append(os.getcwd())
import aero_evo as evo
import objective_funcs as objf


class Gwing_Gui: 

    def __init__(self, root): 

        root.title("G-WING Control Panel")
        root.geometry("700x400")
        self.root = root 
        self.is_running = False 
        self.is_paused = False
        self.seed_wing = None
        settings = json.load(open("tkinter/user_settings.json", "r"))

        # Creating Menubar
        self.menubar = Menu(root)
        file = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label ='File', menu=file)
        file.add_command(label="New Seed", command=self.open_seed_builder)
        file.add_command(label = 'Load Seed', command=None)
        file.add_command(label ='Save Wing', command=None)
        file.add_command(label = "Export Wing", command=None)
        file.add_separator()
        file.add_command(label ='Exit', command=root.destroy)
        edit = Menu(self.menubar, tearoff = 0)
        self.menubar.add_cascade(label ='Edit', menu=edit)
        edit.add_command(label ='Optimizer Settings', command=self.open_optimizer_settings)
        edit.add_command(label ='Wing Parameters', command=self.open_wing_parameters)
        edit.add_separator()
        edit.add_command(label ='XFOIL Settings', command=self.open_xfoil_settings)
        edit.add_command(label ='AVL Settings', command=self.open_avl_settings)

        #creating setup area 
        Label(root, text = "Setup").grid(row=0, column=0, columnspan=3)
        self.fitness_text = StringVar()
        self.fitness_text.set(settings["optimizer settings"]["fitness func"])
        Label(root, textvariable=self.fitness_text).grid(row=1, column=1)
        Label(root, text = "Objective Function:").grid(row=1, column=0)
        b_fitfunc = Button(root, text="Change", command=self.open_function_selector)
        b_fitfunc.grid(row=1, column=2)

        #creating output window 
        Label(root, text = "Output").grid(row=4, column=0, columnspan=3)
        self.text_area = scrolledtext.ScrolledText(root, wrap = WORD, height=10)
        self.text_area.insert(END, "Welcome to G-Wing! The genetic wing optimizer.")
        self.text_area.configure(state="disabled")
        self.text_area.grid(row=5, column=0, columnspan=3)

        #optimizer control buttons
        Label(root, text="Controls").grid(row=2, column=0, columnspan=3)
        self.b_start = Button(root, text="Run", command=self.start_optimizer)
        self.b_pause = Button(root, text="Pause", state=DISABLED, command=self.pause_optimizer)
        self.b_stop = Button(root, text="Stop", state=DISABLED, command=self.stop_optimizer)
        [b.grid(row=3, column=i) for i,b in enumerate([self.b_start,self.b_pause,self.b_stop])]

    def write_to_output(self, text):
        """
        writes text to the output window 
        """ 
        self.text_area.configure(state="normal")
        self.text_area.insert(END, "\n" + text)
        self.text_area.see(END)
        self.text_area.configure(state="disabled")

    def open_optimizer_settings(self):

        self.write_to_output("Opening Optimizer Settings")
        Optimizer_Settings(root)

    def open_wing_parameters(self):
        """
        """
        self.write_to_output("Opening Wing Parameters")
        Wing_Parameters(root)

    def open_xfoil_settings(self):
        """
        """

    def open_avl_settings(self):
        """
        """

    def open_seed_builder(self):
        """
        """

    def open_function_selector(self):
        """
        """
        Function_Selector(self)
        self.write_to_output("Opening Objective Functions")

    def load_inputs(self):
        """
        loads up inputs from user settings and stores as attributes
        """
        full_inputs = json.load(open("tkinter/user_settings.json", "r"))
        self.study_parameters = full_inputs["optimizer settings"]
        self.wing_parameters = full_inputs["wing parameters"]

        fitness_translator = {
            "max L/D":          objf.max_LtoD,
            "max L/D at CL":    objf.LtoD_at_CL,
            "max L^(3/2)/D":    objf.max_L32toD,
            "max L^(1/2)/D":    objf.max_L12toD,
            "max CL max":       objf.max_CL_max
        }

        self.fitness_function = fitness_translator[self.study_parameters["fitness func"]] #! placeholder sample fitness function

    def run_optimizer(self):
        """
        """
        live_plot = True  #! placeholder 
        multiproc = True #! placeholder
        self.load_inputs()

        #initialize population:  
        t_start = time.time()/60
        if live_plot: 
            import matplotlib.pyplot as plt
            fig, axs = evo.initialize_plot()

        #intialize population of m individuals
        popSize = self.study_parameters["population size"]
        k = self.study_parameters["children per generation"]
        if k >= popSize: 
            raise ValueError("invalid children per generation and population size")

        p_gene_mut, p_child_mut = self.study_parameters["gene mutation probability"],\
                                    self.study_parameters["child mutation probability"]

        #clear existing xfoil polar files
        file_list = os.listdir(os.getcwd())
        [os.remove(file) for file in file_list if file.split("_")[0] == "polar"]

        #create population 
        self.population = evo.Population(size=popSize, wing_parameters=self.wing_parameters,\
                                mutation_probs=[p_child_mut, p_gene_mut],\
                                    seed_wing=self.seed_wing, multiproc=multiproc)  

        #evaluate fitness of initial population 
        [chrom.evaluate_fitness(self.fitness_function) for chrom in self.population.chroms]
        self.population.get_best_chrom()

        nIter=None
        if "number of gens" in list(self.study_parameters.keys()):
            nIter = self.study_parameters["number of gens"]
        prev_best = self.population.best_chrom

        #main loop:
        n = 0 # iteration counter 
        while self.is_running: 
            if not self.is_paused: 
                
                if live_plot: evo.update_plot(fig, axs, self.population, prev_best, "running...",\
                                   time0=t_start)
                prev_best = self.population.best_chrom
                self.population.selection(k) #kill off worst performing chromosomes
                self.population.crossover_and_mutation(k) #generate new children 
                [chrom.evaluate_fitness(self.fitness_function) for chrom in self.population.chroms] #evaluate fitness 
                self.population.get_best_chrom()
                n += 1  

            if nIter is not None:
                if n > nIter: self.is_running = False 

    def pause_optimizer(self):
        """
        pauses the optimizer loop
        """
        
        if self.is_paused:
            self.write_to_output("Resuming Solution")
            self.is_paused = False
            self.b_pause.config(text="Pause")
        else:
            self.write_to_output("Pausing solution at next iteration.")
            self.is_paused = True
            self.b_pause.config(text="Resume")

    def stop_optimizer(self):
        """
        stops the optimizer loop
        """
        self.write_to_output("Stopping solution at next iteration.")
        self.is_running = False
        self.b_start.config(state=NORMAL)
        self.b_pause.config(state=DISABLED)
        self.b_stop.config(state=DISABLED)

    def start_optimizer(self):
        """
        continues the optimizer after being paused 
        """
        self.write_to_output("Starting solution")
        self.is_running = True
        self.is_paused = False

        self.b_start.config(state=DISABLED)
        self.b_pause.config(state=NORMAL)
        self.b_stop.config(state=NORMAL)

        # Start a separate thread for the loop
        threading.Thread(target=self.run_optimizer).start() 


class Optimizer_Settings:

    def __init__(self,root):
        
        self.full_settings = json.load(open("tkinter/user_settings.json", "r"))
        self.settings = self.full_settings["optimizer settings"]        
        top = Toplevel(root)
        top.geometry("360x160")
        top.title("optimizer settings")

        Label(top, text="population size:").grid(row=1, column=0)
        self.t_popsize = Entry(top)
        self.t_popsize.insert(END, str(self.settings["population size"]))
        self.t_popsize.grid(row=1, column=1, columnspan=2)    

        Label(top, text="# of children:").grid(row=2, column=0)
        self.t_numChild = Entry(top)
        self.t_numChild.insert(END, str(self.settings["children per generation"]))
        self.t_numChild.grid(row=2, column=1, columnspan=2)

        Label(top, text="gene mutation probability").grid(row=3, column=0)
        self.t_gmut = Entry(top)
        self.t_gmut.insert(END, str(self.settings["gene mutation probability"]))
        self.t_gmut.grid(row=3, column=1, columnspan=2)

        Label(top, text="child mutation probability").grid(row=4, column=0)
        self.t_cmut = Entry(top)
        self.t_cmut.insert(END, str(self.settings["child mutation probability"]))
        self.t_cmut.grid(row=4, column=1, columnspan=2)

        Label(top, text="limit number of generations").grid(row=5, column=0)
        self.numgen_enabled = BooleanVar()
        self.numgen_enabled.set("number of gens" in list(self.settings.keys()))
        cb_numgen = Checkbutton(top, variable=self.numgen_enabled, command=self.num_gens_checkbutton)
        cb_numgen.grid(row=5, column=1)
        self.t_numgen = Entry(top)
        if self.numgen_enabled.get(): 
            self.t_numgen.insert(END, str(self.settings["number of gens"]))
        else: 
            self.t_numgen.config(state=DISABLED)
        self.t_numgen.grid(row=5, column=2)

        self.b_save = Button(master=top, text="save",command=self.save_optimizer_settings)
        self.b_save.grid(row=6, column=0)
        self.b_close = Button(master=top, text="close", command=top.destroy)
        self.b_close.grid(row=6, column=1) 

        self.l_status = Label(top, text="")
        self.l_status.grid(row=7, column=0)

    def save_optimizer_settings(self): 
        
        subdict = {}
        try:
            subdict["fitness func"] = self.settings["fitness func"]
            subdict["population size"] = int(self.t_popsize.get())
            subdict["children per generation"] = int(self.t_numChild.get())
            subdict["gene mutation probability"] = float(self.t_gmut.get())
            subdict["child mutation probability"] = float(self.t_cmut.get())

            if self.numgen_enabled.get():
                subdict["number of gens"] = int(self.t_numgen.get())
        except: 
            self.l_status.configure(text="invalid value")
            return 
        
        if subdict == self.settings: 
            self.l_status.configure(text="settings not changed")
            return 

        self.l_status.configure(text="saved settings.")
        self.settings = subdict
        self.full_settings["optimizer settings"] = self.settings
        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(self.full_settings, file)

    def num_gens_checkbutton(self):

        if self.numgen_enabled.get(): 
            self.t_numgen.config(state=NORMAL) 
        else:
            self.t_numgen.delete(0, END)
            self.t_numgen.config(state=DISABLED) 


class Wing_Parameters: 

    def __init__(self, root):
        """
        
        """
        self.full_settings = json.load(open("tkinter/user_settings.json", "r"))
        self.settings = self.full_settings["wing parameters"]
        top = Toplevel(root) 
        top.geometry("700x350")
        top.title("Wing Parameters")

        Label(top, text="Variable").grid(row=0, column=1)
        Label(top, text="Min").grid(row=0, column=2)
        Label(top, text="Value").grid(row=0, column=3)
        Label(top, text="Max").grid(row=0, column=4)


        # REYNOLDS NUMBER INPUT: 
        init_val = self.settings["Re/chord"]
        Label(top, text="Reynolds Number to Chord:").grid(row=1, column=0)
        self.t_Re_min, self.t_Re_val, self.t_Re_max = Entry(top), Entry(top), Entry(top)
        self.Re_sweep = BooleanVar()
        self.Re_sweep.set(isinstance(init_val, list))
        self.cb_Re = Checkbutton(top, variable=self.Re_sweep, 
                                command=lambda:self.check_button(self.Re_sweep, 
                                                                self.t_Re_min, 
                                                                self.t_Re_val, 
                                                                self.t_Re_max))
        self.cb_Re.grid(row=1, column=1)
        if isinstance(init_val, list):
            self.t_Re_val.config(state=DISABLED)
            self.t_Re_min.insert(END, min(init_val))
            self.t_Re_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_Re_min, self.t_Re_max]]
            self.t_Re_val.insert(END, init_val)
        [t.grid(row=1, column=2+i) for i,t in enumerate([self.t_Re_min, self.t_Re_val, self.t_Re_max])]


        # ASPECT RATIO INPUT
        init_val = self.settings["aspect ratio"]
        Label(top, text="Aspect Ratio:").grid(row=2, column=0)
        self.t_AR_min, self.t_AR_val, self.t_AR_max = Entry(top), Entry(top), Entry(top)
        self.AR_sweep = BooleanVar()
        self.AR_sweep.set(isinstance(init_val, list))
        self.cb_AR = Checkbutton(top, variable=self.AR_sweep, 
                                command=lambda:self.check_button(self.AR_sweep, 
                                                                self.t_AR_min, 
                                                                self.t_AR_val, 
                                                                self.t_AR_max))
        self.cb_AR.grid(row=2, column=1)
        if isinstance(init_val, list):
            self.t_AR_val.config(state=DISABLED)
            self.t_AR_min.insert(END, min(init_val))
            self.t_AR_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_AR_min, self.t_AR_max]]
            self.t_AR_val.insert(END, init_val)
        [t.grid(row=2, column=2+i) for i,t in enumerate([self.t_AR_min, self.t_AR_val, self.t_AR_max])]  


        # TAPER RATIO INPUT 
        init_val = self.settings["taper"]
        Label(top, text="Taper Ratio:").grid(row=3, column=0)
        self.t_tap_min, self.t_tap_val, self.t_tap_max = Entry(top), Entry(top), Entry(top)
        self.tap_sweep = BooleanVar()
        self.tap_sweep.set(isinstance(init_val, list))
        self.cb_tap = Checkbutton(top, variable=self.tap_sweep, 
                                command=lambda:self.check_button(self.tap_sweep, 
                                                                self.t_tap_min, 
                                                                self.t_tap_val, 
                                                                self.t_tap_max))
        self.cb_tap.grid(row=3, column=1)
        if isinstance(init_val, list):
            self.t_tap_val.config(state=DISABLED)
            self.t_tap_min.insert(END, min(init_val))
            self.t_tap_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_tap_min, self.t_tap_max]]
            self.t_tap_val.insert(END, init_val)
        [t.grid(row=3, column=2+i) for i,t in enumerate([self.t_tap_min, self.t_tap_val, self.t_tap_max])]  


        # QUARTER CHORD SWEEP INPUT 
        init_val = self.settings["sweep deg"]
        Label(top, text="Q-Chord Sweep (deg):").grid(row=4, column=0)
        self.t_sweep_min, self.t_sweep_val, self.t_sweep_max = Entry(top), Entry(top), Entry(top)
        self.sweep_sweep = BooleanVar()
        self.sweep_sweep.set(isinstance(init_val, list))
        self.cb_sweep = Checkbutton(top, variable=self.sweep_sweep, 
                                    command=lambda:self.check_button(self.sweep_sweep, 
                                                                    self.t_sweep_min, 
                                                                    self.t_sweep_val, 
                                                                    self.t_sweep_max))
        self.cb_sweep.grid(row=4, column=1)
        if isinstance(init_val, list):
            self.t_sweep_val.config(state=DISABLED)
            self.t_sweep_min.insert(END, min(init_val))
            self.t_sweep_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_sweep_min, self.t_sweep_max]]
            self.t_sweep_val.insert(END, init_val)
        [t.grid(row=4, column=2+i) for i,t in enumerate([self.t_sweep_min, self.t_sweep_val, self.t_sweep_max])] 


        # TIP TWIST INPUT 
        init_val = self.settings["twist deg"]
        Label(top, text="Tip Twist (deg):").grid(row=5, column=0)
        self.t_twist_min, self.t_twist_val, self.t_twist_max = Entry(top), Entry(top), Entry(top)
        self.twist_sweep = BooleanVar()
        self.twist_sweep.set(isinstance(init_val, list))
        self.cb_twist = Checkbutton(top, variable=self.twist_sweep, 
                                    command=lambda:self.check_button(self.twist_sweep, 
                                                                    self.t_twist_min, 
                                                                    self.t_twist_val,
                                                                    self.t_twist_max))
        self.cb_twist.grid(row=5, column=1)
        if isinstance(init_val, list):
            self.t_twist_val.config(state=DISABLED)
            self.t_twist_min.insert(END, min(init_val))
            self.t_twist_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_twist_min, self.t_twist_max]]
            self.t_twist_val.insert(END, init_val)
        [t.grid(row=5, column=2+i) for i,t in enumerate([self.t_twist_min, self.t_twist_val, self.t_twist_max])]


        Label(top, text="Root Airfoil:").grid(row=6, column=0, columnspan=2)
        # ROOT AIRFOIL CAMBER 
        init_val = self.settings["root camber"]
        Label(top, text="1st NACA digit:").grid(row=7, column=0)
        self.t_rootm_min, self.t_rootm_val, self.t_rootm_max = Entry(top), Entry(top), Entry(top)
        self.rootm_sweep = BooleanVar()
        self.rootm_sweep.set(isinstance(init_val, list))
        self.cb_rootm = Checkbutton(top, variable=self.rootm_sweep, 
                                    command=lambda:self.check_button(self.rootm_sweep,
                                                                     self.t_rootm_min,
                                                                     self.t_rootm_val,
                                                                     self.t_rootm_max))
        self.cb_rootm.grid(row=7, column=1)
        if isinstance(init_val, list):
            self.t_rootm_val.config(state=DISABLED)
            self.t_rootm_min.insert(END, min(init_val))
            self.t_rootm_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_rootm_min, self.t_rootm_max]]
            self.t_rootm_val.insert(END, init_val)
        [t.grid(row=7, column=2+i) for i,t in enumerate([self.t_rootm_min, self.t_rootm_val, self.t_rootm_max])]


        # ROOT AIRFOIL CAMBER LOCATION
        init_val = self.settings["root camber location"]
        Label(top, text="2nd NACA digit:").grid(row=8, column=0)
        self.t_rootp_min, self.t_rootp_val, self.t_rootp_max = Entry(top), Entry(top), Entry(top)
        self.rootp_sweep = BooleanVar()
        self.rootp_sweep.set(isinstance(init_val, list))
        self.cb_rootp = Checkbutton(top, variable=self.rootp_sweep,
                                    command=lambda:self.check_button(self.rootp_sweep,
                                                                     self.t_rootp_min,
                                                                     self.t_rootp_val,
                                                                     self.t_rootp_max))
        self.cb_rootp.grid(row=8, column=1)
        if isinstance(init_val, list):
            self.t_rootp_val.config(state=DISABLED)
            self.t_rootp_min.insert(END, min(init_val))
            self.t_rootp_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_rootp_min, self.t_rootp_max]]
            self.t_rootp_val.insert(END, init_val)
        [t.grid(row=8, column=2+i) for i,t in enumerate([self.t_rootp_min, self.t_rootp_val, self.t_rootp_max])]


        # ROOT AIRFOIL THICKNESS
        init_val = self.settings["root camber location"]
        Label(top, text="3+4th NACA digits:").grid(row=9, column=0)
        self.t_roott_min, self.t_roott_val, self.t_roott_max = Entry(top), Entry(top), Entry(top)
        self.roott_sweep = BooleanVar()
        self.roott_sweep.set(isinstance(init_val, list))
        self.cb_roott = Checkbutton(top, variable=self.roott_sweep,
                                    command=lambda:self.check_button(self.roott_sweep,
                                                                     self.t_roott_min,
                                                                     self.t_roott_val,
                                                                     self.t_roott_max))
        self.cb_roott.grid(row=9, column=1)
        if isinstance(init_val, list):
            self.t_roott_val.config(state=DISABLED)
            self.t_roott_min.insert(END, min(init_val))
            self.t_roott_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_roott_min, self.t_roott_max]]
            self.t_roott_val.insert(END, init_val)
        [t.grid(row=9, column=2+i) for i,t in enumerate([self.t_roott_min, self.t_roott_val, self.t_roott_max])]

        Label(top, text="Tip Airfoil:").grid(row=10, column=0, columnspan=2)
        
        
        # ROOT AIRFOIL CAMBER 
        init_val = self.settings["tip camber"]
        Label(top, text="1st NACA digit:").grid(row=11, column=0)
        self.t_tipm_min, self.t_tipm_val, self.t_tipm_max = Entry(top), Entry(top), Entry(top)
        self.tipm_sweep = BooleanVar()
        self.tipm_sweep.set(isinstance(init_val, list))
        self.cb_tipm = Checkbutton(top, variable=self.tipm_sweep,
                                   command=lambda:self.check_button(self.tipm_sweep,
                                                                    self.t_tipm_min,
                                                                    self.t_tipm_val,
                                                                    self.t_tipm_max))
        self.cb_tipm.grid(row=11, column=1)
        if isinstance(init_val, list):
            self.t_tipm_val.config(state=DISABLED)
            self.t_tipm_min.insert(END, min(init_val))
            self.t_tipm_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_tipm_min, self.t_tipm_max]]
            self.t_tipm_val.insert(END, init_val)
        [t.grid(row=11, column=2+i) for i,t in enumerate([self.t_tipm_min, self.t_tipm_val, self.t_tipm_max])]


        # ROOT AIRFOIL CAMBER LOCATION
        init_val = self.settings["tip camber location"]
        Label(top, text="2nd NACA digit:").grid(row=12, column=0)
        self.t_tipp_min, self.t_tipp_val, self.t_tipp_max = Entry(top), Entry(top), Entry(top)
        self.tipp_sweep = BooleanVar()
        self.tipp_sweep.set(isinstance(init_val, list))
        self.cb_tipp = Checkbutton(top, variable=self.tipp_sweep,
                                   command=lambda:self.check_button(self.tipp_sweep,
                                                                    self.t_tipp_min,
                                                                    self.t_tipp_val,
                                                                    self.t_tipp_max))
        self.cb_tipp.grid(row=12, column=1)
        if isinstance(init_val, list):
            self.t_tipp_val.config(state=DISABLED)
            self.t_tipp_min.insert(END, min(init_val))
            self.t_tipp_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_tipp_min, self.t_tipp_max]]
            self.t_tipp_val.insert(END, init_val)
        [t.grid(row=12, column=2+i) for i,t in enumerate([self.t_tipp_min, self.t_tipp_val, self.t_tipp_max])]


        # ROOT AIRFOIL THICKNESS
        init_val = self.settings["tip camber location"]
        Label(top, text="3+4th NACA digits:").grid(row=13, column=0)
        self.t_tipt_min, self.t_tipt_val, self.t_tipt_max = Entry(top), Entry(top), Entry(top)
        self.tipt_sweep = BooleanVar()
        self.tipt_sweep.set(isinstance(init_val, list))
        self.cb_tipt = Checkbutton(top, variable=self.tipt_sweep,
                                   command=lambda:self.check_button(self.tipt_sweep,
                                                                    self.t_tipt_min,
                                                                    self.t_tipt_val,
                                                                    self.t_tipt_max))
        self.cb_tipt.grid(row=13, column=1)
        if isinstance(init_val, list):
            self.t_tipt_val.config(state=DISABLED)
            self.t_tipt_min.insert(END, min(init_val))
            self.t_tipt_max.insert(END, max(init_val))
        elif isinstance(init_val, (int, float)):
            [t.config(state=DISABLED) for t in [self.t_tipt_min, self.t_tipt_max]]
            self.t_tipt_val.insert(END, init_val)
        [t.grid(row=13, column=2+i) for i,t in enumerate([self.t_tipt_min, self.t_tipt_val, self.t_tipt_max])]        


        Button(master=top, text="save", command=self.save_wing_parameters).grid(row=14, column=0)
        Button(master=top, text="close", command=top.destroy).grid(row=14, column=1)
        self.l_status = Label(top, text="")
        self.l_status.grid(row=15, column=0)

    def save_wing_parameters(self):
        """
        
        """
        subdict = {}
        def set_dict_entry(key, boolvar, min_field, val_field, max_field):
            if boolvar.get():
                subdict[key] = [float(min_field.get()), float(max_field.get())]
            else: 
                subdict[key] = float(val_field.get())

        try: 
            set_dict_entry("aspect ratio", self.AR_sweep, self.t_AR_min, self.t_AR_val, self.t_AR_max)
            set_dict_entry("taper", self.tap_sweep, self.t_tap_min, self.t_tap_val, self.t_tap_max)
            set_dict_entry("sweep deg", self.sweep_sweep, self.t_sweep_min, self.t_sweep_val, self.t_sweep_max)
            set_dict_entry("twist deg", self.twist_sweep, self.t_twist_min, self.t_twist_val, self.t_twist_max)
            set_dict_entry("root camber", self.rootm_sweep, self.t_rootm_min, self.t_rootm_val, self.t_rootm_max)
            set_dict_entry("root camber location", self.rootp_sweep, self.t_rootp_min, self.t_rootp_val, self.t_rootp_max)
            set_dict_entry("root thickness", self.roott_sweep, self.t_roott_min, self.t_roott_val, self.t_roott_max)
            set_dict_entry("tip camber", self.tipm_sweep, self.t_tipm_min, self.t_tipm_val, self.t_tipm_max)
            set_dict_entry("tip camber location", self.tipp_sweep, self.t_tipp_min, self.t_tipp_val, self.t_tipp_max)
            set_dict_entry("tip thickness", self.tipt_sweep, self.t_tipt_min, self.t_tipt_val, self.t_tipt_max)
            set_dict_entry("Re/chord", self.Re_sweep, self.t_Re_min, self.t_Re_val, self.t_Re_max)
        except: 
            self.l_status.configure(text="invalid value")
            return

        #check if subdict has been changed
        if self.settings == subdict: 
            self.l_status.configure(text="settings not changed")
            return 

        self.l_status.configure(text="saved settings")
        self.full_settings["wing parameters"] = subdict
        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(self.full_settings, file)

    def check_button(self, var, min_field, val_field, max_field):
        """
        """
        if var.get(): #is list 
            min_field.config(state=NORMAL)
            val_field.config(state=DISABLED)
            max_field.config(state=NORMAL)
        else: 
            min_field.config(state=DISABLED)
            val_field.config(state=NORMAL)
            max_field.config(state=DISABLED)


class Function_Selector:

    def __init__(self, mainGui):
        
        self.main = mainGui 
        top = Toplevel(mainGui.root) 
        top.geometry("300x300")
        top.title("Objective Functions")
        self.items=["max L/D", "max L/D at CL","max L^(3/2)/D","max L^(1/2)/D","max CL max"]

        self.selected_item = None
        self.listbox = Listbox(top, selectmode=SINGLE)
        self.listbox.grid(row=0, column=0, columnspan=2)
        for item in self.items:
            self.listbox.insert(END, item)

        self.b_select = Button(top, text="select", command=self.select_item)
        self.b_select.grid(row=2, column=0)

        self.b_close = Button(top, text="close", command=lambda:top.destroy())
        self.b_close.grid(row=2, column=1)

        self.b_custom = Button(top, text="use custom", command=None)
        self.b_custom.grid(row=1, column=0, columnspan=2)

        self.l_status = Label(top, text="")
        self.l_status.grid(row=3, column=0)

    def select_item(self):

        selected_index = self.listbox.curselection()
        if selected_index:
            self.selected_item = self.items[selected_index[0]]

        full_settings = json.load(open("tkinter/user_settings.json", "r"))
        full_settings["optimizer settings"]["fitness func"] = self.selected_item

        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(full_settings, file)

        self.l_status.configure(text="selected")
        self.main.fitness_text.set(self.selected_index)
        self.main.root.update()
            

if __name__=="__main__":
    
    root = Tk()
    gui = Gwing_Gui(root)
    root.config(menu=gui.menubar)
    root.mainloop()
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
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import time


class Gwing_Gui: 

    def __init__(self, root): 

        root.title("G-WING Control Panel")
        root.geometry("1200x900")
        self.root = root 
        self.is_running = False 
        self.is_paused = False
        self.seed_wing = None
        self.settings = json.load(open("tkinter/user_settings.json", "r"))

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

        #creating tabs 
        self.tabs = Notebook(root)
        self.tab_main = Frame(self.tabs)
        self.tab_geom = Frame(self.tabs)
        self.tab_xfoil = Frame(self.tabs)
        self.tab_avl = Frame(self.tabs)
        self.tabs.add(self.tab_main, text="Main")
        self.tabs.add(self.tab_geom, text="Geometry")
        self.tabs.add(self.tab_xfoil, text="XFOIL")
        self.tabs.add(self.tab_avl, text="AVL")
        self.tabs.grid(row=0, column=0, columnspan=3)

        #initialize figures and tabs
        self.initialize_figures()
        self.initialize_geom_tab()
        self.initialize_main_tab()
        self.initialize_xfoil_tab()
        self.initialize_avl_tab()

        #creating output window 
        Label(root, text = "Output").grid(row=4, column=0, columnspan=3)
        self.text_area = scrolledtext.ScrolledText(root, wrap = WORD, height=10)
        self.text_area.insert(END, "Welcome to G-Wing! The genetic wing optimizer.")
        self.text_area.configure(state="disabled")
        self.text_area.grid(row=5, column=0, columnspan=3)

        #optimizer control buttons
        Label(root, text="Controls").grid(row=2, column=0, columnspan=3)
        self.frame_controls = Frame(root)
        self.b_start = Button(self.frame_controls, text="Run", command=self.start_optimizer)
        self.b_pause = Button(self.frame_controls, text="Pause", state=DISABLED, command=self.pause_optimizer)
        self.b_stop = Button(self.frame_controls, text="Stop", state=DISABLED, command=self.stop_optimizer)
        [b.grid(row=0, column=i) for i,b in enumerate([self.b_start,self.b_pause,self.b_stop])]
        self.frame_controls.grid(row=3, column=0, columnspan=3)

    def initialize_main_tab(self):
        """
        """
        canvas = FigureCanvasTkAgg(self.fig_hist, master=self.tab_main)
        canvas.draw() 
        canvas.get_tk_widget().grid(row=2, column=0, columnspan=3, rowspan=2)

        #creating setup area 
        self.status = StringVar()
        self.status.set("waiting")
        self.frame_setup = Frame(self.tab_main)
        Label(self.frame_setup, text = "Setup").grid(row=0, column=0)
        Label(self.frame_setup, text="Status: ").grid(row=0, column=3)
        Label(self.frame_setup, textvariable=self.status).grid(row=0, column=4)
        Label(self.frame_setup, text="Runtime: ").grid(row=1, column=3)
        self.timer_label = Label(self.frame_setup, text="00:00:00")
        self.timer_label.grid(row=1, column=4)
        self.fitness_text = StringVar()
        self.fitness_text.set(self.settings["optimizer settings"]["fitness func"])
        Label(self.frame_setup, textvariable=self.fitness_text).grid(row=1, column=1)
        Label(self.frame_setup, text = "Objective Function:  ").grid(row=1, column=0)
        b_fitfunc = Button(self.frame_setup, text="Change", command=self.open_function_selector)
        b_fitfunc.grid(row=1, column=2)
        self.frame_setup.grid(row=0, column=0, columnspan=3)

    def initialize_geom_tab(self):
        """
        """
        canvas = FigureCanvasTkAgg(self.fig_geom, master=self.tab_geom)
        canvas.draw() 
        canvas.get_tk_widget().grid(row=0, column=0, columnspan=3)

    def initialize_xfoil_tab(self):
        """
        """
        root = self.tab_xfoil
        canvas = FigureCanvasTkAgg(self.fig_xfoil, master=root)
        canvas.draw()
        canvas.get_tk_widget().grid(row=2, column=0, columnspan=10)

        self.b_runxfoil = Button(root, text="Run XFOIL", command=None).grid(row=0, column=0)
        Label(root, text="Airfoil 1: NACA").grid(row=0, column=1)
        Label(root, text="Airfoil 2 NACA").grid(row=1, column=1)
        self.t_foil1 = Entry(root).grid(row=0, column=2)
        self.t_foil2 = Entry(root).grid(row=1, column=2)

        Label(root, text="XFOIL Parameters").grid(row=4, column=0)

        Label(root, text="Iterations:").grid(row=5, column=0)
        self.t_niter = Entry(root).grid(row=5, column=1)

        Label(root, text="Ncrit:").grid(row=5, column=2)
        self.t_ncrit = Entry(root).grid(row=5, column=3)

        Label(root, text="Alpha Sequence:").grid(row=5, column=4)
        Label(root, text="initial").grid(row=4, column=5)
        Label(root, text="step").grid(row=4, column=6)
        Label(root, text="final").grid(row=4, column=7)
        self.t_alpha_i = Entry(root).grid(row=5, column=5)
        self.t_alpha_step = Entry(root).grid(row=5, column=6)
        self.t_alpha_f = Entry(root).grid(row=5, column=7)

        Label(root, text="Reynolds/Chord:").grid(row=5, column=8)
        self.t_re_c = Entry(root).grid(row=5, column=9)


    def initialize_avl_tab(self):
        """
        """
        root = self.tab_avl

        root = self.tab_avl
        canvas = FigureCanvasTkAgg(self.fig_avl, master=root)
        canvas.draw()
        canvas.get_tk_widget().grid(row=2, column=0, columnspan=7)

        self.frame_avl_plots = Frame(root)
        self.b_runavl = Button(root, text="Run AVL", command=None).grid(row=0, column=0)

        Label(root, text="Aspect Ratio:").grid(row=0, column=1)
        self.t_AR = Entry(root).grid(row=0, column=2)
        Label(root, text="Taper Ratio:").grid(row=1, column=1)
        self.t_taper = Entry(root).grid(row=1, column=2)
        Label(root, text="Sweep (deg):").grid(row=0, column=3)
        self.t_sweep = Entry(root).grid(row=0, column=4)
        Label(root, text="Twist (deg):").grid(row=1, column=3)
        self.t_twist = Entry(root).grid(row=1, column=4)
        Label(root, text="Root: NACA").grid(row=0, column=5)
        self.t_rootfoil = Entry(root).grid(row=0, column=6)
        Label(root, text="Tip: NACA").grid(row=1, column=5)
        self.t_tipfoil = Entry(root).grid(row=1, column=6)

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
        self.write_to_output("Opening xfoil settings")
        Xfoil_Settings(root)

    def open_avl_settings(self):
        """
        """
        self.write_to_output("Opening avl settings")
        Avl_Settings(root)

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
        self.xfoil_settings = full_inputs["xfoil settings"]
        self.avl_settings = full_inputs["avl settings"]

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
        Runs the optimizer population initialization and iterative loop
        """
    
        self.load_inputs()
        multiproc = self.study_parameters["num cpu"]
        
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
        self.write_to_output(text="Initializing population...")
        self.population = evo.Population(size=popSize, 
                            wing_parameters=self.wing_parameters,
                            mutation_probs=[p_child_mut, p_gene_mut],
                            xfoil_settings=self.xfoil_settings, 
                            avl_settings=self.avl_settings, 
                            seed_wing=None, multiproc=multiproc)  

        #evaluate fitness of initial population 
        [chrom.evaluate_fitness(self.fitness_function) for chrom in self.population.chroms]
        self.population.get_best_chrom()

        nIter=None
        if "number of gens" in list(self.study_parameters.keys()):
            nIter = self.study_parameters["number of gens"]
        prev_best = self.population.best_chrom

        #main loop:
        n = 0 # iteration counter 
        self.write_to_output(text="Iterating population...")
        while self.is_running: 
            if not self.is_paused: 
                
                self.update_history_plot(self.population)
                self.update_geometry_plots(self.population, prev_best)
                #evo.update_3figplot(self.fig_live, self.axs_live, self.population, prev_best, "running...",\
                                   #time0=t_start)
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
            threading.Thread(target=self.update_timer).start()
        else:
            self.write_to_output("Pausing solution at next iteration.")
            self.is_paused = True
            self.b_pause.config(text="Resume")
        
        self.status.set("paused")

    def stop_optimizer(self):
        """
        stops the optimizer loop
        """
        self.write_to_output("Stopping solution at next iteration.")

        [ax.clear() for ax in [self.ax_foils, self.ax_planf, self.ax_hist]]

        self.is_running = False
        self.b_start.config(state=NORMAL)
        self.b_pause.config(state=DISABLED)
        self.b_stop.config(state=DISABLED)
        self.status.set("stopped")

    def start_optimizer(self):
        """
        continues the optimizer after being paused 
        """
        self.write_to_output("Starting solution")
        self.is_running = True
        self.is_paused = False
        threading.Thread(target=self.start_timer).start()
        self.status.set("running")

        self.b_start.config(state=DISABLED)
        self.b_pause.config(state=NORMAL)
        self.b_stop.config(state=NORMAL)

        # Start a separate thread for the loop
        threading.Thread(target=self.run_optimizer).start() 

    def initialize_figures(self):
        """
        initializes history and geometry plots 
        """
        #history plot 
        self.fig_hist = plt.figure(figsize=(12,4), dpi=100)
        self.ax_hist = self.fig_hist.add_subplot(111)
        self.ax_hist.spines[:].set_color("dimgrey")
        self.ax_hist.set_xlabel("Iteration")
        self.ax_hist.set_ylabel("Objective")
        self.ax_hist.grid(color="lightgrey")
    
        #planform and airfoils plot
        self.fig_geom = plt.figure(figsize=(12,6), dpi=100)
        self.ax_planf = self.fig_geom.add_subplot(211)
        self.ax_planf.set_title("Planform", fontdict={"color":"grey"})
        self.ax_foils = self.fig_geom.add_subplot(212)
        self.ax_foils.set_title("Root & Tip Sections", fontdict={"color":"grey"})
        self.ax_foils.spines[:].set_color("lightgrey")
        self.ax_foils.tick_params(colors="grey")
        self.ax_planf.spines[:].set_color("lightgrey")
        self.ax_planf.tick_params(colors="grey")

        #xfoil plots 
        self.fig_xfoil = plt.figure(figsize=(12,5), dpi=100)
        self.ax_xfoil_cla = self.fig_xfoil.add_subplot(131)
        self.ax_xfoil_cma = self.fig_xfoil.add_subplot(132)
        self.ax_xfoil_cdcl = self.fig_xfoil.add_subplot(133)

        #avl plots
        self.fig_avl = plt.figure(figsize=(12,5), dpi=100)

    def update_history_plot(self, population):
        """
        updates the chromosome fitness history plot 
        """
        if hasattr(self, "line_hist"):
            self.line_hist.remove()
        n_list = list(range(len(population.best_chrom_fitness)))

        self.line_hist, = self.ax_hist.plot(n_list,population.best_chrom_fitness, "-o", color="red",\
             linewidth=0.5, markersize=3)
        
        ver_count = 1
        self.ax_hist.annotate(text="init", xy=(0, population.best_chrom_fitness[0]), ha="right", va="bottom")
        for i,f in enumerate(population.best_chrom_fitness[:-1]):
            if population.best_chrom_fitness[i+1] > f: 
                xy = (i+1, population.best_chrom_fitness[i+1])
                self.ax_hist.annotate(text=str(ver_count), xy=xy, ha="right", va="bottom")
                ver_count+=1

        self.fig_hist.canvas.draw()
        self.fig_hist.canvas.flush_events()
    
    def update_geometry_plots(self, population, prev_best):
        """
        updates planform and geometry plots
        """
        
        if hasattr(self, "lines_planf"):
            [l.remove() for l in self.lines_planf]

        if hasattr(self, "lines_foils"):
            [l.remove() for l in self.lines_foils]

        best_wing = population.best_chrom
        lines1 = prev_best.plot_wing_planform(self.ax_planf, linecolor="lightgrey")
        lines2 = best_wing.plot_wing_planform(self.ax_planf, linecolor="black")
        lines3 = prev_best.plot_wing_airfoils(self.ax_foils, linecolor="lightgrey")
        lines4 = best_wing.plot_wing_airfoils(self.ax_foils, linecolor="black")

        self.lines_planf = lines1 + lines2
        self.lines_foils = lines3 + lines4

        self.fig_geom.canvas.draw() 
        self.fig_geom.canvas.flush_events()

    def start_timer(self):
        global start_time
        start_time = time.time() 
        self.update_timer()

    def update_timer(self):
        current_time = time.time() - start_time
        hours = int(current_time // 3600)
        minutes = int((current_time // 60) % 60)
        seconds = int(current_time % 60)
        self.timer_label.config(text=f"{hours:02d}:{minutes:02d}:{seconds:02d}")
        if not self.is_paused or not self.is_running:
            self.timer_label.after(1000, self.update_timer)

    def reset_timer(self):
        self.timer_label.config(text="00:00:00")


class Optimizer_Settings:

    def __init__(self,root):
        
        self.full_settings = json.load(open("tkinter/user_settings.json", "r"))
        self.settings = self.full_settings["optimizer settings"]        
        top = Toplevel(root)
        top.geometry("360x180")
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

        Label(top, text="# of parallel processors").grid(row=5, column=0)
        self.t_numcor = Entry(top)
        self.t_numcor.insert(END, str(self.settings["num cpu"]))
        self.t_numcor.grid(row=5, column=1, columnspan=2)

        Label(top, text="limit number of generations").grid(row=6, column=0)
        self.numgen_enabled = BooleanVar()
        self.numgen_enabled.set("number of gens" in list(self.settings.keys()))
        cb_numgen = Checkbutton(top, variable=self.numgen_enabled, command=self.num_gens_checkbutton)
        cb_numgen.grid(row=6, column=1)
        self.t_numgen = Entry(top)
        if self.numgen_enabled.get(): 
            self.t_numgen.insert(END, str(self.settings["number of gens"]))
        else: 
            self.t_numgen.config(state=DISABLED)
        self.t_numgen.grid(row=6, column=2)

        self.b_save = Button(master=top, text="save",command=self.save_optimizer_settings)
        self.b_save.grid(row=7, column=0)
        self.b_close = Button(master=top, text="close", command=top.destroy)
        self.b_close.grid(row=7, column=1) 

        self.l_status = Label(top, text="")
        self.l_status.grid(row=8, column=0)

    def save_optimizer_settings(self): 
        
        subdict = {}
        try:
            subdict["fitness func"] = self.settings["fitness func"]
            subdict["population size"] = int(self.t_popsize.get())
            subdict["children per generation"] = int(self.t_numChild.get())
            subdict["gene mutation probability"] = float(self.t_gmut.get())
            subdict["child mutation probability"] = float(self.t_cmut.get())
            subdict["num cpu"] = int(self.t_numcor.get())

            if self.numgen_enabled.get():
                subdict["number of gens"] = int(self.t_numgen.get())
        except: 
            self.l_status.configure(text="invalid value")
            return 
        
        if subdict["num cpu"] < 1: 
            self.l_status.configure(text="invalid number of cores")
            return
        
        if subdict == self.settings: 
            self.l_status.configure(text="settings not changed")
            return 

        self.l_status.configure(text="saved settings")
        self.settings = subdict
        self.full_settings["optimizer settings"] = self.settings
        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(self.full_settings, file)

    def num_gens_checkbutton(self):
        #number of generations check button controller
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
        init_val = self.settings["root thickness"]
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
        init_val = self.settings["tip thickness"]
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
        self.main.fitness_text.set(self.selected_item)
        self.main.root.update()
            

class Xfoil_Settings: 

    def __init__(self, root):

        self.full_settings = json.load(open("tkinter/user_settings.json", "r"))
        self.settings = self.full_settings["xfoil settings"]
        top = Toplevel(root) 
        top.geometry("360x180")
        top.title("Default XFOIL Settings")

        Label(top, text="Iterations:").grid(row=1, column=0)
        self.t_iter = Entry(top)
        self.t_iter.insert(END, str(self.settings["iterations"]))
        self.t_iter.grid(row=1, column=1, columnspan=2) 

        Label(top, text="Transition Parameter (ncrit):").grid(row=2, column=0)
        self.t_ncrit = Entry(top)
        self.t_ncrit.insert(END, str(self.settings["Ncrit"]))
        self.t_ncrit.grid(row=2, column=1, columnspan=2) 

        Label(top, text="Timeout Limit (sec):").grid(row=3, column=0)
        self.t_timeout = Entry(top)
        self.t_timeout.insert(END, str(self.settings["timeout limit"]))
        self.t_timeout.grid(row=3, column=1, columnspan=2) 

        Label(top, text="Initial AoA (deg):").grid(row=4, column=0)
        self.t_alpha_i = Entry(top)
        self.t_alpha_i.insert(END, str(self.settings["alpha_i"]))
        self.t_alpha_i.grid(row=4, column=1, columnspan=2) 

        Label(top, text="Final AoA (deg):").grid(row=5, column=0)
        self.t_alpha_f = Entry(top)
        self.t_alpha_f.insert(END, str(self.settings["alpha_f"]))
        self.t_alpha_f.grid(row=5, column=1, columnspan=2) 

        Label(top, text="AoA Step (deg):").grid(row=6, column=0)
        self.t_alpha_step = Entry(top)
        self.t_alpha_step.insert(END, str(self.settings["alpha_step"]))
        self.t_alpha_step.grid(row=6, column=1, columnspan=2) 

        self.b_save = Button(master=top, text="save",command=self.save_xfoil_settings)
        self.b_save.grid(row=7, column=0)
        self.b_close = Button(master=top, text="close", command=top.destroy)
        self.b_close.grid(row=7, column=1) 

        self.l_status = Label(top, text="")
        self.l_status.grid(row=8, column=0)

    def save_xfoil_settings(self):

        subdict = {}
        try:
            subdict["iterations"] = int(self.t_iter.get())
            subdict["Ncrit"] = int(self.t_ncrit.get())
            subdict["timeout limit"] = float(self.t_timeout.get())
            subdict["alpha_i"] = float(self.t_alpha_i.get())
            subdict["alpha_f"] = float(self.t_alpha_f.get())
            subdict["alpha_step"] = float(self.t_alpha_step.get())

        except: 
            self.l_status.configure(text="invalid value")
            return 
        
        if subdict == self.settings: 
            self.l_status.configure(text="settings not changed")
            return 

        self.l_status.configure(text="saved settings")
        self.settings = subdict
        self.full_settings["xfoil settings"] = self.settings
        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(self.full_settings, file)


class Avl_Settings: 

    def __init__(self, root):

        self.full_settings = json.load(open("tkinter/user_settings.json", "r"))
        self.settings = self.full_settings["avl settings"]
        top = Toplevel(root) 
        top.geometry("360x200")
        top.title("Default AVL Settings")

        Label(top, text="# Chord Vortices:").grid(row=1, column=0)
        self.t_nchord = Entry(top)
        self.t_nchord.insert(END, str(self.settings["Nchord"]))
        self.t_nchord.grid(row=1, column=1, columnspan=2) 

        Label(top, text="# Halfspan Vortices:").grid(row=2, column=0)
        self.t_nspan = Entry(top)
        self.t_nspan.insert(END, str(self.settings["Nspan"]))
        self.t_nspan.grid(row=2, column=1, columnspan=2) 

        Label(top, text="Span Spacing:").grid(row=3, column=0)
        self.t_sspace = Entry(top)
        self.t_sspace.insert(END, str(self.settings["Sspace"]))
        self.t_sspace.grid(row=3, column=1, columnspan=2) 

        Label(top, text="Chord Spacing:").grid(row=4, column=0)
        self.t_cspace = Entry(top)
        self.t_cspace.insert(END, str(self.settings["Cspace"]))
        self.t_cspace.grid(row=4, column=1, columnspan=2) 

        Label(top, text="Initial AoA (deg):").grid(row=5, column=0)
        self.t_alpha_i = Entry(top)
        self.t_alpha_i.insert(END, str(self.settings["alpha_i"]))
        self.t_alpha_i.grid(row=5, column=1, columnspan=2) 

        Label(top, text="Final AoA (deg):").grid(row=6, column=0)
        self.t_alpha_f = Entry(top)
        self.t_alpha_f.insert(END, str(self.settings["alpha_f"]))
        self.t_alpha_f.grid(row=6, column=1, columnspan=2) 

        Label(top, text="AoA Step (deg):").grid(row=7, column=0)
        self.t_alpha_step = Entry(top)
        self.t_alpha_step.insert(END, str(self.settings["alpha_step"]))
        self.t_alpha_step.grid(row=7, column=1, columnspan=2) 

        self.b_save = Button(master=top, text="save",command=self.save_xfoil_settings)
        self.b_save.grid(row=8, column=0)
        self.b_close = Button(master=top, text="close", command=top.destroy)
        self.b_close.grid(row=8, column=1) 

        self.l_status = Label(top, text="")
        self.l_status.grid(row=9, column=0)

    def save_xfoil_settings(self):

        subdict = {}
        try:
            subdict["Nchord"] = int(self.t_nchord.get())
            subdict["Nspan"] = int(self.t_nspan.get())
            subdict["Sspace"] = int(self.t_sspace.get())
            subdict["Cspace"] = int(self.t_cspace.get())
            subdict["alpha_i"] = float(self.t_alpha_i.get())
            subdict["alpha_f"] = float(self.t_alpha_f.get())
            subdict["alpha_step"] = float(self.t_alpha_step.get())

        except: 
            self.l_status.configure(text="invalid value")
            return

        if subdict["Sspace"] not in [-3,-2,-1,0,1,2,3]:
            self.l_status.configure(text="invalid span spacing")
            return 

        if subdict["Cspace"] not in [-3,-2,-1,0,1,2,3]:
            self.l_status.configure(text="invalid chord spacing")
            return 
        
        if subdict == self.settings: 
            self.l_status.configure(text="settings not changed")
            return 

        self.l_status.configure(text="saved settings")
        self.settings = subdict
        self.full_settings["avl settings"] = self.settings
        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(self.full_settings, file)


if __name__ == "__main__":
    
    root = Tk()
    gui = Gwing_Gui(root)
    root.config(menu=gui.menubar)
    root.mainloop()
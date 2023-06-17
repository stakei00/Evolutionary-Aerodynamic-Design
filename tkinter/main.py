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
        self.write_to_output("Opening Objective Functions")
        Function_Selector()

    def load_inputs(self):
        """
        loads up inputs from user settings and stores as attributes
        """
        full_inputs = json.load(open("tkinter/user_settings.json", "r"))
        self.study_parameters = full_inputs["optimizer settings"]
        self.wing_parameters = full_inputs["wing parameters"]
        self.fitness_function = objf.max_LtoD #! placeholder sample fitness function

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

    def save_optimizer_settings(self): 
        
        subdict = {}
        subdict["population size"] = int(self.t_popsize.get())
        subdict["children per generation"] = int(self.t_numChild.get())
        subdict["gene mutation probability"] = float(self.t_gmut.get())
        subdict["child mutation probability"] = float(self.t_cmut.get())

        if self.numgen_enabled.get():
            subdict["number of gens"] = int(self.t_numgen.get())

        if subdict == self.settings: 
            print("settings not changed")
            return 

        
        print("saved settings.")
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
        full_settings = json.load(open("tkinter/user_settings.json", "r"))
        settings = full_settings["wing parameters"]
        top = Toplevel(root) 
        top.geometry("700x350")
        top.title("Wing Parameters")

        Label(top, text="Variable").grid(row=0, column=1)
        Label(top, text="Min").grid(row=0, column=2)
        Label(top, text="Value").grid(row=0, column=3)
        Label(top, text="Max").grid(row=0, column=4)

        # REYNOLDS NUMBER INPUT: 
        init_val = settings["Re/chord"]
        Label(top, text="Reynolds Number to Chord:").grid(row=1, column=0)
        self.Re_sweep = BooleanVar()
        self.Re_sweep.set(isinstance(init_val, list))
        self.cb_Re = Checkbutton(top, variable=self.Re_sweep)
        self.cb_Re.grid(row=1, column=1)
        self.t_Re_min, self.t_Re_val, self.t_Re_max = Entry(top), Entry(top), Entry(top)
        if isinstance(init_val, list):
            self.t_Re_min.insert(END, min(init_val))
            self.t_Re_max.insert(END, max(init_val))
        elif isinstance(init_val, int or float):
            self.t_Re_val.insert(END, init_val)
        [t.grid(row=1, column=2+i) for i,t in enumerate([self.t_Re_min, self.t_Re_val, self.t_Re_max])]

        # ASPECT RATIO INPUT
        init_val = settings["aspect ratio"]
        Label(top, text="Aspect Ratio:").grid(row=2, column=0)
        self.AR_sweep = BooleanVar()
        self.AR_sweep.set(isinstance(init_val, list))
        self.cb_AR = Checkbutton(top, variable=self.AR_sweep)
        self.cb_AR.grid(row=2, column=1)
        self.t_AR_min, self.t_AR_val, self.t_AR_max = Entry(top), Entry(top), Entry(top)
        if isinstance(init_val, list):
            self.t_AR_min.insert(END, min(init_val))
            self.t_AR_max.insert(END, max(init_val))
        elif isinstance(init_val, int or float):
            self.t_AR_val.insert(END, init_val)
        [t.grid(row=2, column=2+i) for i,t in enumerate([self.t_AR_min, self.t_AR_val, self.t_AR_max])]  

        # TAPER RATIO INPUT 
        init_val = settings["taper"]
        Label(top, text="Taper Ratio:").grid(row=3, column=0)
        self.tap_sweep = BooleanVar()
        self.tap_sweep.set(isinstance(init_val, list))
        self.cb_tap = Checkbutton(top, variable=self.tap_sweep)
        self.cb_tap.grid(row=3, column=1)
        self.t_tap_min, self.t_tap_val, self.t_tap_max = Entry(top), Entry(top), Entry(top)
        if isinstance(init_val, list):
            self.t_tap_min.insert(END, min(init_val))
            self.t_tap_max.insert(END, max(init_val))
        elif isinstance(init_val, int or float):
            self.t_tap_val.insert(END, init_val)
        [t.grid(row=3, column=2+i) for i,t in enumerate([self.t_tap_min, self.t_tap_val, self.t_tap_max])]  


        # QUARTER CHORD SWEEP INPUT 
        init_val = settings["sweep deg"]
        Label(top, text="Q-Chord Sweep (deg):").grid(row=4, column=0)
        self.sweep_sweep = BooleanVar()
        self.sweep_sweep.set(isinstance(init_val, list))
        self.cb_sweep = Checkbutton(top, variable=self.sweep_sweep)
        self.cb_sweep.grid(row=4, column=1)
        self.t_sweep_min, self.t_sweep_val, self.t_sweep_max = Entry(top), Entry(top), Entry(top)
        if isinstance(init_val, list):
            self.t_sweep_min.insert(END, min(init_val))
            self.t_sweep_max.insert(END, max(init_val))
        elif isinstance(init_val, int or float):
            self.t_sweep_val.insert(END, init_val)
        [t.grid(row=4, column=2+i) for i,t in enumerate([self.t_sweep_min, self.t_sweep_val, self.t_sweep_max])] 

        Label(top, text="Tip Twist (deg)").grid(row=5, column=0)
        t_twist = Entry(top)
        t_twist.insert(END, str(settings["twist deg"]))
        t_twist.grid(row=5, column=1)

        Label(top, text="Root airfoil:").grid(row=6, column=0, columnspan=2)
        Label(top, text="1st NACA digit:").grid(row=7, column=0)
        t_root_m = Entry(top)
        t_root_m.insert(END, str(settings["root camber"]))
        t_root_m.grid(row=7, column=1)

        Label(top, text="2nd NACA digit:").grid(row=8, column=0)
        t_root_p = Entry(top)
        t_root_p.insert(END, str(settings["root camber location"]))
        t_root_p.grid(row=8, column=1)

        Label(top, text="3+4th NACA digits").grid(row=9, column=0)
        t_root_t = Entry(top)
        t_root_t.insert(END, str(settings["root camber"]))
        t_root_t.grid(row=9, column=1)

        Label(top, text="Tip airfoil:").grid(row=10, column=0, columnspan=2)
        Label(top, text="1st NACA digit:").grid(row=11, column=0)
        t_tip_m = Entry(top)
        t_tip_m.insert(END, str(settings["tip camber"]))
        t_tip_m.grid(row=11, column=1)

        Label(top, text="2nd NACA digit:").grid(row=12, column=0)
        t_tip_p = Entry(top)
        t_tip_p.insert(END, str(settings["tip camber location"]))
        t_tip_p.grid(row=12, column=1)

        Label(top, text="3+4th NACA digits").grid(row=13, column=0)
        t_tip_t = Entry(top)
        t_tip_t.insert(END, str(settings["tip camber"]))
        t_tip_t.grid(row=13, column=1)

        Button(master=top, text="save", 
               command=lambda:(self.save_wing_parameters(full_settings), 
                self.write_to_output("saved wing parameters"))).grid(row=14, column=0)

        Button(master=top, text="close", command=top.destroy).grid(row=14, column=1)

    def save_wing_parameters(self):
        """
        
        """
        return #! debugging make it work for lists, not just ints/floats 
        subdict = {}
        subdict["aspect ratio"] = float(t_AR.get())
        subdict["taper"] = float(t_taper.get())
        subdict["sweep deg"] = float(t_sweep.get())
        subdict["twist deg"] = float(t_twist.get())
        subdict["root camber"] = float(t_root_m.get())/100
        subdict["root camber location"] = float(t_root_p.get())/10
        subdict["root thickness"] = float(t_root_t.get())/100
        subdict["tip camber"] = float(t_tip_m.get())/100
        subdict["tip camber location"] = float(t_tip_p.get())/10
        subdict["tip thickness"] = float(t_tip_t.get())/100
        subdict["Re/chord"] = round(int(t_Re_c),-5)
        settings["wing parameters"] = subdict
        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(settings, file)


class Function_Selector:

    def __init__(self):
        top = Toplevel(root) 
        top.geometry("700x350")
        top.title("Objective Functions")

if __name__=="__main__":
    
    root = Tk()
    gui = Gwing_Gui(root)
    root.config(menu=gui.menubar)
    root.mainloop()
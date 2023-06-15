import sys, os
sys.path.append(os.getcwd())
import aero_evo as evo 
from tkinter import *
from tkinter.ttk import *
from tkinter import scrolledtext
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from gwing_commands import*
import json 


#main parameters 
root = Tk()
root.title("G-WING")
root.geometry("1200x900")

#child window functions: 
def open_optimizer_settings():

    def save_optimizer_settings(settings):
    
        subdict = {}
        subdict["population size"] = int(t_popsize.get())
        subdict["children per generation"] = int(t_numChild.get())
        subdict["gene mutation probability"] = float(t_gmut.get())
        subdict["child mutation probability"] = float(t_cmut.get())
        subdict["number of gens"] = int(t_numgen.get())

        settings["Optimizer Settings"] = subdict

        with open("tkinter/user_settings.json", "w") as file: 
            json.dump(settings, file)

    full_settings = json.load(open("tkinter/user_settings.json", "r"))
    settings = full_settings["optimizer settings"]
    top = Toplevel(root) 
    top.geometry("270x170")
    top.title("optimizer settings")
    
    Label(top, text="population size:").grid(row=1, column=0)
    t_popsize = Entry(top)
    t_popsize.insert(END, str(settings["population size"]))
    t_popsize.grid(row=1, column=1)    

    Label(top, text="# of children:").grid(row=2, column=0)
    t_numChild = Entry(top)
    t_numChild.insert(END, str(settings["children per generation"]))
    t_numChild.grid(row=2, column=1)

    Label(top, text="gene mutation probability").grid(row=3, column=0)
    t_gmut = Entry(top)
    t_gmut.insert(END, str(settings["gene mutation probability"]))
    t_gmut.grid(row=3, column=1)

    Label(top, text="child mutation probability").grid(row=4, column=0)
    t_cmut = Entry(top)
    t_cmut.insert(END, str(settings["child mutation probability"]))
    t_cmut.grid(row=4, column=1)

    Label(top, text="number of generations").grid(row=5, column=0)
    t_numgen = Entry(top)
    t_numgen.insert(END, str(settings["number of gens"]))
    t_numgen.grid(row=5, column=1)

    Button(master=top, text="save", 
           command=lambda:(save_optimizer_settings(full_settings), 
            write_to_output("saved optimizer settings"))).grid(row=6, column=0)
    
    Button(master=top, text="close", command=top.destroy).grid(row=6, column=1)

def open_wing_parameters(): 

    def save_wing_parameters(settings):
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

    full_settings = json.load(open("tkinter/user_settings.json", "r"))
    settings = full_settings["wing parameters"]
    top = Toplevel(root) 
    top.geometry("300x300")
    top.title("Wing Parameters")

    Label(top, text="Reynolds Number to Chord:").grid(row=0, column=0)
    t_Re_c = Entry(top)
    t_Re_c.insert(END, str(settings["Re/chord"]))
    t_Re_c.grid(row=0, column=1)  
    
    Label(top, text="Aspect Ratio:").grid(row=1, column=0)
    t_AR = Entry(top)
    t_AR.insert(END, str(settings["aspect ratio"]))
    t_AR.grid(row=1, column=1)    

    Label(top, text="Taper Ratio:").grid(row=2, column=0)
    t_taper = Entry(top)
    t_taper.insert(END, str(settings["taper"]))
    t_taper.grid(row=2, column=1)

    Label(top, text="Q-Chord Sweep (deg)").grid(row=3, column=0)
    t_sweep = Entry(top)
    t_sweep.insert(END, str(settings["sweep deg"]))
    t_sweep.grid(row=3, column=1)

    Label(top, text="Tip Twist (deg)").grid(row=4, column=0)
    t_twist = Entry(top)
    t_twist.insert(END, str(settings["twist deg"]))
    t_twist.grid(row=4, column=1)

    Label(top, text="Root airfoil:").grid(row=5, column=0, columnspan=2)
    Label(top, text="1st NACA digit:").grid(row=6, column=0)
    t_root_m = Entry(top)
    t_root_m.insert(END, str(settings["root camber"]))
    t_root_m.grid(row=6, column=1)

    Label(top, text="2nd NACA digit:").grid(row=7, column=0)
    t_root_p = Entry(top)
    t_root_p.insert(END, str(settings["root camber location"]))
    t_root_p.grid(row=7, column=1)

    Label(top, text="3+4th NACA digits").grid(row=8, column=0)
    t_root_t = Entry(top)
    t_root_t.insert(END, str(settings["root camber"]))
    t_root_t.grid(row=8, column=1)

    Label(top, text="Tip airfoil:").grid(row=9, column=0, columnspan=2)
    Label(top, text="1st NACA digit:").grid(row=10, column=0)
    t_tip_m = Entry(top)
    t_tip_m.insert(END, str(settings["tip camber"]))
    t_tip_m.grid(row=10, column=1)

    Label(top, text="2nd NACA digit:").grid(row=11, column=0)
    t_tip_p = Entry(top)
    t_tip_p.insert(END, str(settings["tip camber location"]))
    t_tip_p.grid(row=11, column=1)

    Label(top, text="3+4th NACA digits").grid(row=12, column=0)
    t_tip_t = Entry(top)
    t_tip_t.insert(END, str(settings["tip camber"]))
    t_tip_t.grid(row=12, column=1)

    Button(master=top, text="save", 
           command=lambda:(save_wing_parameters(full_settings), 
            write_to_output("saved wing parameters"))).grid(row=13, column=0)
    
    Button(master=top, text="close", command=top.destroy).grid(row=13, column=1)

def open_xfoil_settings():

    settings = json.load(open("tkinter/user_settings.json", "r"))["xfoil settings"]
    top = Toplevel(root) 
    top.geometry("500x500")
    top.title("xfoil settings")
    Label(top, text="edit xfoil settings here").grid(row=0, column=0)

def open_avl_settings():

    settings = json.load(open("tkinter/user_settings.json", "r"))["avl settings"]
    top = Toplevel(root) 
    top.geometry("500x500")
    top.title("avl")
    Label(top, text="edit AVL settings here").grid(row=0, column=0) 

# Creating Menubar
menubar = Menu(root)
file = Menu(menubar, tearoff = 0)
menubar.add_cascade(label ='File', menu = file)
file.add_command(label="New Seed", command=None)
file.add_command(label = 'Load Seed', command = None)
file.add_command(label ='Save Wing', command = None)
file.add_command(label = "Export Wing", command=None)
file.add_separator()
file.add_command(label ='Exit', command = root.destroy)
edit = Menu(menubar, tearoff = 0)
menubar.add_cascade(label ='Edit', menu = edit)
edit.add_command(label ='Optimizer Settings', command=open_optimizer_settings)
edit.add_command(label ='Wing Parameters', command=open_wing_parameters)
edit.add_separator()
edit.add_command(label ='XFOIL Settings', command=open_xfoil_settings)
edit.add_command(label ='AVL Settings', command=open_avl_settings)

# Adding tabs to notebook
tabControl = Notebook(root)
tab1 = Frame(tabControl)
tab2 = Frame(tabControl)
tab3 = Frame(tabControl)
tabControl.add(tab1, text ='Main')
tabControl.add(tab2, text ='Geometry')
tabControl.add(tab3, text ='Aerodynamics')
tabControl.grid(row=0, column=0, columnspan=5 ,sticky="nsew")

#adding output to bottom 
def write_to_output(text):
    text_area.configure(state="normal")
    text_area.insert(END, "\n" + text)
    text_area.see(END)
    text_area.configure(state="disabled")

Label(root, text = "Output").grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
text_area = scrolledtext.ScrolledText(root, wrap = WORD, height=10)
text_area.insert(END, "Welcome to G-Wing! The genetic wing optimizer")
text_area.configure(state="disabled")
text_area.grid(row=2, column=0)

#optimizer control buttons
b_start = Button(root, text="start new", command=lambda:(write_to_output("new solution started"), 
                                                    b_start.configure(state="disabled"),
                                                    b_stop.configure(state="normal"),
                                                    b_pause.configure(state="normal"),
                                                    run_optimizer()))

b_pause = Button(root, text="pause", command=lambda: (write_to_output("solution paused"),
                                                      b_play.configure(state="normal"),
                                                      b_pause.configure(state="disabled")))

b_play = Button(root, text="play", command=lambda:(write_to_output("continuing solution"),
                                                   b_play.configure(state="disabled"),
                                                   b_pause.configure(state="normal")))

b_stop = Button(root, text="stop", command=lambda:(write_to_output("solution stopped"),
                                                   b_start.configure(state="normal"),
                                                   b_pause.configure(state="disabled"),
                                                   b_play.configure(state="disabled"),
                                                   b_stop.configure(state="disabled")))
#initial button states:
b_play.configure(state="disabled")
b_pause.configure(state="disabled")
b_stop.configure(state="disabled")
[b.grid(row=2, column=1+i) for i,b in enumerate([b_start,b_pause,b_play,b_stop])]

#main tab: 
summ_frame = Frame(tab1, width=800, height=200)
hist_frame = Frame(tab1, width=1200, height=400)
Label(tab1, text="Summary").grid(row=0, column=0)
Label(tab1, text="Genetic History").grid(row=2, column=0)
summ_frame.grid(row=1, column=0)
hist_frame.grid(row=3, column=0)

fig,axs = evo.initialize_plot()
canvas = FigureCanvasTkAgg(fig, master=tab1)
canvas.draw()
canvas.get_tk_widget().grid(row=3,column=0)

#geometry tab 
plan_frame = Frame(tab2, width=1200, height=200)
sect_frame = Frame(tab2, width=1200, height=200)
Label(tab2, text="Planform").grid(row=0, column=0)
Label(tab2, text="Sections").grid(row=2, column=0)
plan_frame.grid(row=1, column=0)
sect_frame.grid(row=3, column=0)

if __name__ == "__main__":
    root.config(menu=menubar)
    mainloop()
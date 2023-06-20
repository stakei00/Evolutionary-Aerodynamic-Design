import math 
import numpy as np 

class Wing: 

    def __init__(self, airfoil_root: object, airfoil_tip: object, \
                taper_ratio:float, aspect_ratio:float, sweep_deg: float, \
                    twist_deg:float) -> None: 
        """
        """
        self.airfoil_root = airfoil_root
        self.airfoil_tip = airfoil_tip
        self.Re_c = airfoil_root.Re #reynolds # to root chord ratio 
        self.tip_twist_deg = twist_deg
        self.sweep_deg = sweep_deg
        self.AR = aspect_ratio
        self.taper = taper_ratio
        self.root_chord = 1
        self.tip_chord = self.taper*self.root_chord
        self.span = self.AR*(self.root_chord+self.tip_chord)/2
        self.area = self.span**2/self.AR
        self.mean_aero_chord = 2*self.root_chord*((self.taper**2 + self.taper +\
                                                    1)/(self.taper + 1))/3

        self.run_avl_on_wing(-6, 20, 1)
        self.interpolate_avl_results()
        self.get_max_lift_coefficient()

    def run_avl_on_wing(self, alpha_i:float, alpha_f:float,\
                        alpha_step:float) -> None: 
        """
        runs avl on wing configuration
        """

        import avlpy.avlInput as avlInput
        
        #setting up header
        header = avlInput.AvlHeader(configname="WING", iYsym=1, Sref=self.area,\
                                    Cref=self.mean_aero_chord, Bref=self.span)
        
        #setting up root section 
        section_root = avlInput.AvlSection(chord=self.root_chord, Nspan=10, \
                            afile=["NACA",self.airfoil_root.NACA_4series_desig],\
                                CDCL=self.airfoil_root.CDCL_avl,Sspace=1)

        #setting up tip section
        self.x_t = self.root_chord/4 - self.tip_chord/4 + \
            math.tan(math.radians(self.sweep_deg))*self.span/2
        
        section_tip = avlInput.AvlSection(le=(self.x_t, self.span/2, 0), \
                            chord=self.tip_chord, ainc=self.tip_twist_deg, \
                                Nspan = 10, afile = \
                                    ["NACA",self.airfoil_tip.NACA_4series_desig],\
                                        CDCL=self.airfoil_tip.CDCL_avl, Sspace=1)

        #setting up main wing surface 
        surface = avlInput.AvlSurface(name="Wing", Nchord=12, Nspan=20, Sspace=1, Cspace=1)
        surface.addSections([section_root, section_tip])

        #running 
        avlInputObj = avlInput.AvlInput(header, [surface])
        self.avl_results,_ = avlInput.runAvl(avlInputObj, \
                                alpha=np.arange(alpha_i, alpha_f, alpha_step), \
                                    collect_surface_forces=False)

    def interpolate_avl_results(self) -> None:

        self.alpha = np.linspace(self.avl_results.Alpha[0], \
                                            self.avl_results.Alpha[-1], 100)
        self.lift_coefficient = np.interp(self.alpha, self.avl_results.Alpha, \
                                            self.avl_results.CLtot)
        self.drag_coefficient = np.interp(self.alpha, self.avl_results.Alpha, \
                                            self.avl_results.CDtot) 
        self.moment_coefficient = np.interp(self.alpha, self.avl_results.Alpha,\
                                            self.avl_results.Cmtot)     

    def get_max_lift_coefficient(self) -> None: 
        """
        analyzes lift coefficient distribution along wing and determines 
        maximum lift coefficient based on root and tip airfoil
        self.run_avl_on_wing() must be called beforehand 
        """
        self.max_lift_coefficient = None
        alpha = self.avl_results.Alpha
        lift_coefficient = self.avl_results.CLtot
        y_le = self.avl_results.stripForces[0].strips[0].yle

        spanwise_cl = []
        max_spanwise_cl = []
        max_span_cl_y = []
        for stripForces in self.avl_results.stripForces: 
            spanwise_cl.append(stripForces.strips[0].cl)
            i_max_cl  = np.argmax(stripForces.strips[0].cl) #index of max cl
            max_spanwise_cl.append(stripForces.strips[0].cl[i_max_cl])
            max_span_cl_y.append(y_le[i_max_cl])

        for i,a in enumerate(alpha):
            max_cl = max_spanwise_cl[i]
            y_max_cl = max_span_cl_y[i]
            cl_max_root = self.airfoil_root.max_lift_coefficient
            cl_max_tip = self.airfoil_tip.max_lift_coefficient
            m = (cl_max_tip - cl_max_root)/(self.span/2)

            limit_max_cl = m*y_max_cl + cl_max_root

            if limit_max_cl <= max_cl:
                self.max_lift_coefficient = lift_coefficient[i]
                break 

    def plot_wing_planform(self, ax, linecolor:str="grey") -> None:
        """
        
        """
        horz_ax_data = [0,0,0, self.span/2,self.span/2, self.span/2,self.span/2, 0]
        vert_ax_data = [0, -self.root_chord,-self.root_chord, \
                        -self.x_t-self.tip_chord,-self.x_t-self.tip_chord, \
                            -self.x_t,-self.x_t, 0]
        
        line1, = ax.plot(horz_ax_data, vert_ax_data, color=linecolor, linewidth=1)
        line2, = ax.plot([-1*x for x in horz_ax_data], vert_ax_data, color=linecolor, linewidth=1)
        line3, = ax.plot([-self.span/2, 0, self.span/2],[-self.x_t-self.tip_chord/4,\
                        -self.root_chord/4, -self.x_t-self.tip_chord/4],\
                            color=linecolor, linewidth=0.75, linestyle="dashed")
        #ax.fill(horz_ax_data, vert_ax_data, facecolor="aliceblue")
        #ax.fill([-1*x for x in horz_ax_data], vert_ax_data, facecolor="aliceblue")

        ax.set_aspect("equal", adjustable="box")
        return [line1, line2, line3]

    def plot_wing_airfoils(self, ax, linecolor:str="grey") -> None:
        """
        
        """
        line1 = self.airfoil_root.plot_airfoil(ax, chord_scale=self.root_chord, \
                                       linecolor=linecolor)
        line2 = self.airfoil_tip.plot_airfoil(ax, chord_scale=self.tip_chord, \
                                twist_deg=self.tip_twist_deg, \
                                    origin=(-self.x_t,0), linecolor=linecolor)
        
        ax.set_aspect("equal", adjustable="box")
        return [line1, line2]
    
    def plot_avl_xfoil_results(self) -> None: 
        """
        creates a figure with several subplots showing results from AVL and xfoil 
        """
        import matplotlib.pyplot as plt 
        fig = plt.figure(figsize=(12,8))
        fig.canvas.manager.set_window_title("AVL XFOIL Results")
        fig.suptitle("AVL and XFOIL Results")
        gs = fig.add_gridspec(3,3)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[0,2])
        ax4 = fig.add_subplot(gs[1,0])
        ax5 = fig.add_subplot(gs[1,1])
        ax6 = fig.add_subplot(gs[2,:2])
        ax7 = fig.add_subplot(gs[1:,2])

        #airfoil plots: 
        ax1.plot(self.airfoil_root.alpha, self.airfoil_root.lift_coefficient,\
                 label=f"r: NACA{self.airfoil_root.NACA_4series_desig}")
        ax1.plot(self.airfoil_tip.alpha, self.airfoil_tip.lift_coefficient,\
                 label=f"t: NACA{self.airfoil_tip.NACA_4series_desig}")
        ax1.set_xlabel("$\u03B1^{\circ}$"), ax1.set_ylabel("$c_l$"), ax1.grid()
        ax1.legend()

        ax2.plot(self.airfoil_root.alpha, self.airfoil_root.moment_coefficient)
        ax2.plot(self.airfoil_tip.alpha, self.airfoil_tip.moment_coefficient)
        ax2.set_xlabel("$\u03B1^{\circ}$"), ax2.set_ylabel("$c_m$"), ax2.grid()

        ax3.plot(self.airfoil_root.drag_coefficient, self.airfoil_root.lift_coefficient)
        ax3.plot(self.airfoil_tip.drag_coefficient, self.airfoil_tip.lift_coefficient)
        ax3.set_xlabel("$c_d$"), ax3.set_ylabel("$c_l$"), ax3.grid()
        
        #wing plots: 
        ax4.plot(self.alpha, self.lift_coefficient, color="red", label="wing")
        ax4.set_xlabel("$\u03B1^{\circ}$"), ax4.set_ylabel("$C_L$"), ax4.grid()
        ax4.legend()

        ax5.plot(self.alpha, self.moment_coefficient, color="red")
        ax5.set_xlabel("$\u03B1^{\circ}$"), ax5.set_ylabel("$C_M$"), ax5.grid()

        ax7.plot(self.drag_coefficient, self.lift_coefficient, color="red")
        ax7.set_xlabel("$C_D$"), ax7.set_ylabel("$C_L$"), ax7.grid()

        for i,a in enumerate(self.avl_results.Alpha):
            if a % 2 != 0: continue
            cl = self.avl_results.stripForces[i].strips[0].cl
            y = self.avl_results.stripForces[i].strips[0].yle
            ax6.plot(y, cl, label=f"\u03B1 = {a}", linewidth=0.75, color="red")
            ax6.annotate(str(a), (y[0], cl[0]), ha="right", va="center")
        
        ax6.plot([0, self.span/2],[self.airfoil_root.max_lift_coefficient,\
                                    self.airfoil_tip.max_lift_coefficient],\
                                        linestyle="dashed", color="black")
        ax6.set_xlabel("y"), ax6.set_ylabel("$c_l$"), ax6.grid()
        

        plt.show()

    def export_wing(self): 
        import pickle
        import datetime

        now = datetime.now()
        current_time = now.strftime("%H_%M")
        with open(f"wing_{current_time}", "wb") as file:
            pickle.dump(self, file)
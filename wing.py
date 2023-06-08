import math 
import numpy as np 

class Wing: 

    def __init__(self, airfoil_root: object, airfoil_tip: object, span: float, \
                taper_ratio:float, aspect_ratio:float, sweep_deg: float, \
                    twist_deg:float) -> None: 
        """
        """
        self.airfoil_root = airfoil_root
        self.airfoil_tip = airfoil_tip
        self.tip_twist_deg = twist_deg
        self.sweep_deg = sweep_deg
        self.AR = aspect_ratio
        self.taper = taper_ratio
        self.span = span

        self.area = self.span**2/self.AR
        self.root_chord = 2*self.area/(self.span*(1 + self.taper))
        self.tip_chord = self.taper*self.root_chord
        self.mean_aero_chord = 2*self.root_chord*((self.taper**2 + self.taper +\
                                                    1)/(self.taper + 1))/3

        self.run_avl_on_wing(-6, 20, 1)
        self.interpolate_avl_results()
        self.get_max_lift_coefficient()

    def run_avl_on_wing(self, alpha_i:float, alpha_f:float, alpha_step:float) -> None: 
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
                                CDCL=self.airfoil_root.CDCL_avl)

        #setting up tip section
        self.x_t = self.root_chord/4 - self.tip_chord/4 + \
            math.tan(math.radians(self.sweep_deg))*self.span/2
        
        section_tip = avlInput.AvlSection(le=(self.x_t, self.span/2, 0), \
                            chord=self.tip_chord, ainc=self.tip_twist_deg, \
                                Nspan = 10, afile = \
                                    ["NACA",self.airfoil_tip.NACA_4series_desig],\
                                        CDCL=self.airfoil_tip.CDCL_avl)

        #setting up main wing surface 
        surface = avlInput.AvlSurface(name="Wing", Nchord=12, Nspan=20)
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

    def plot_wing_planform(self, ax, linecolor:str="white") -> None:
        """
        
        """
        horz_ax_data = [0,0,0, self.span/2,self.span/2, self.span/2,self.span/2, 0]
        vert_ax_data = [0, -self.root_chord,-self.root_chord, \
                        -self.x_t-self.tip_chord,-self.x_t-self.tip_chord, \
                            -self.x_t,-self.x_t, 0]
        ax.plot(horz_ax_data, vert_ax_data, color=linecolor)
        ax.plot([-1*x for x in horz_ax_data], vert_ax_data, color=linecolor)
        #ax.fill(horz_ax_data, vert_ax_data, facecolor="aliceblue")
        #ax.fill([-1*x for x in horz_ax_data], vert_ax_data, facecolor="aliceblue")
        ax.set_aspect("equal", adjustable="box")

    def plot_wing_airfoils(self, ax, linecolor:str="white") -> None:
        """
        
        """
        self.airfoil_root.plot_airfoil(ax, chord_scale=self.root_chord, \
                                       linecolor=linecolor)
        self.airfoil_tip.plot_airfoil(ax, chord_scale=self.tip_chord, \
                                twist_deg=self.tip_twist_deg, \
                                    origin=(-self.x_t,0), linecolor=linecolor)
        
        ax.set_aspect("equal", adjustable="box")

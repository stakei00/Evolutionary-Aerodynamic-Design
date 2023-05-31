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

        self.run_avl_on_wing()

    def run_avl_on_wing(self) -> None: 
        """
        runs avl on wing configuration
        """

        import avlpy.avlInput as avlInput
        import numpy as np 
        import math 
        
        #setting up header
        header = avlInput.AvlHeader(configname="WING", iYsym=1, Sref=self.area,\
                                    Cref=self.mean_aero_chord, Bref=self.span)
        
        #setting up root section 
        section_root = avlInput.AvlSection(chord=self.root_chord, Nspan=10, \
                            afile="NACA"+self.airfoil_root.NACA_4series_desig)

        #setting up tip section
        x_t = self.root_chord/4 - self.tip_chord/4 + \
            math.tan(math.radians(self.sweep_deg))*self.span/2
        
        section_tip = avlInput.AvlSection(le=(-x_t, self.span/2, 0), \
                            chord=self.tip_chord, ainc=self.tip_twist_deg, \
                                Nspan = 10, afile = \
                                    "NACA"+self.airfoil_tip.NACA_4series_desig)

        #setting up main wing surface 
        surface = avlInput.AvlSurface(name="Wing", Nchord=12, Nspan=20)
        surface.addSections(section_root, section_tip)

        #running 
        avlInputObj = avlInput.AvlInput(header, [surface])
        results,_ = avlInput.runAvl(avlInputObj, alpha=np.linspace(-10, 10, 1))

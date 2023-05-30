import numpy as np 
import math 

class Airfoil:

    def __init__(self, m:float, p:float, t:float, c:float) -> None:
        """
        Initializes airfoil object 
        """
        self.max_camber = m
        self.max_camber_loc = p
        self.thickness_to_chord = t
        self.chord = c
        
        self.NACA_4series_desig = str(int(1000*round(m*100, 0) + 100*round(p*10, 0) + round(t*100, 0)))

    """
    def get_surface_coords(self, m:float, p:float, t:float, nPts:int) -> None:
        
        #generates discrete airfoil surface using cambered naca 4 series equations
        #Inputs: 
        #    m: maximum camber 
        #    p: location of max camber
        #    t: thickness to chord ratio 
        #    nPts: number of points on upper/lower surface (linspace'd) 
        #Example: NACA 2412 uses 
        #    m = 0.02
        #    p = 0.4
        #    t = 0.12
        
        def y_c(x):
            if 0 <= x <= p: 
                return m*(2*p*x - x**2)/(p**2)
            elif p < x <= 1: 
                return m*((1-2*p) + 2*p*x - x**2)/((1-p)**2)
            else: return None 
            
        def dydx_c(x):
            if 0 <= x <= p:
                return 2*m*(p-x)/p**2
            elif p < x <= 1: 
                return 2*m*(p-x)/((1-p)**2)
            else: return None

        def y_t(x):
            return 5*t*(0.2969*math.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**2)


        xList = np.linspace(0,1, nPts)
        
        x_upper, y_upper, x_lower, y_lower = [],[],[],[]
        
        for x in xList:
            thet = math.atan(dydx_c(x))
            
            x_upper.append(x - y_t(x)*math.sin(thet))
            y_upper.append(y_c(x) + y_t(x)*math.cos(thet))

            x_lower.append(x + y_t(x)*math.sin(thet))
            y_lower.append(y_c(x) - y_t(x)*math.cos(thet))

        x_upper.reverse()
        y_upper.reverse()

        self.x = np.array(x_upper + x_lower[1:])
        self.y = np.array(y_upper + y_lower[1:])
    """

    def run_XFOIL(self, alpha_i:float, alpha_f:float, alpha_step:float, \
                  Re:int or float, n_iter:int) -> None: 
        """
        Analysis of airfoil in XFOIL. Adds results to object attributes
        """
        import os
        import subprocess
        import numpy as np

        airfoil_name_naca = "6415"

        if os.path.exists("polar_file.txt"):
            os.remove("polar_file.txt")

        input_file = open("input_file.in", 'w')
        #input_file.write(f"LOAD {airfoil_name}.dat\n")
        input_file.write(f"NACA {airfoil_name_naca}\n")
        #input_file.write(airfoil_name + '\n')
        input_file.write("PANE\n")
        input_file.write("OPER\n")
        input_file.write(f"Visc {Re}\n")
        input_file.write("PACC\n")
        input_file.write("polar_file.txt\n\n")
        input_file.write(f"ITER {n_iter}\n")
        input_file.write(f"ASeq {alpha_i} {alpha_f} {alpha_step}\n")
        input_file.write("\n\n")
        input_file.write("quit\n")
        input_file.close()
        subprocess.call("xfoil.exe < input_file.in", shell=True)

        polar_data = np.loadtxt("polar_file.txt", skiprows=12)

        self.alpha = polar_data[:,0]
        self.lift_coefficient = polar_data[:,1]
        self.drag_coefficient = polar_data[:,2]
        self.moment_coefficient = polar_data[:,4] 

    def find_3pt_drag_polar(self) -> list: 
        """
        Gets 3 points on CLCD drag polar for use with AVL 
        """

        self.CDCL_avl = []

        perc_change_thresh = 0.25
        abs_perc_change = [None]

        for i,Cd in enumerate(self.drag_coefficient[1:]):
            perc_change = (Cd - self.drag_coefficient[i-1])/(0.5*(Cd + \
                                                self.drag_coefficient[i-1]))
            
            abs_perc_change.append(abs(perc_change))

        #find first CLCD point (neg CL region)
        CL_CD = [self.lift_coefficient[0], self.drag_coefficient[0]]
        
        for i, perc in enumerate(abs_perc_change[1:]):
            if self.lift_coefficient[i] > 0: continue
            if perc >= perc_change_thresh: 
                CL_CD[0] = self.lift_coefficient[i]
                CL_CD[1] = self.drag_coefficient[i]
        
        self.CDCL_avl.append(CL_CD[0])
        self.CDCL_avl.append(CL_CD[1])

        #find min CD: 
        idx_CDmin = self.drag_coefficient.argmin()
        self.CDCL_avl.append(self.lift_coefficient[idx_CDmin])
        self.CDCL_avl.append(self.drag_coefficient[idx_CDmin])

        #find last CLCD point (pos CL region)
        for i, perc in enumerate(abs_perc_change[1:]):
            if self.lift_coefficient[i] < 0: continue 

            if perc <= perc_change_thresh:
                CL_CD[0] = self.lift_coefficient[i]
                CL_CD[1] = self.drag_coefficient[i]

        self.CDCL_avl.append(CL_CD[0])
        self.CDCL_avl.append(CL_CD[1])

    def get_max_lift_coefficient(self) -> None: 
        """
        estimates the max lift coefficient of the airfoil 
        """

    
    def set_position_and_angle(self, halfspan:float, x_offset:float, \
                               y_offset:float, inclination_deg:float) -> None: 
        """
        sets the position and angle of the airfoil section. For wing, use RHS 
        coordinate system with roll(x)-axis point out of the nose, pitch(y)-axis
        pointed out of right wing, and yaw(z) axis pointing down. Origin at 
        leading vertex of root airfoil. 
        """
        self.y_wing = halfspan
        self.x_wing = x_offset
        self.z_wing = y_offset
        self.inc_deg = inclination_deg

if __name__ == "__main__":
    m = 0.02
    p = 0.4
    t = 0.12
    c = 1

    af = Airfoil(m,p,t,c)
    af.run_XFOIL(alpha_i=-10, alpha_f=15, alpha_step=1, Re=1e6, n_iter=100)
    af.find_3pt_drag_polar()

    import matplotlib.pyplot as plt 
    plt.figure()
    plt.plot(af.alpha, af.lift_coefficient)
    plt.grid() 
    plt.xlabel("alpha (deg)"), plt.ylabel("C_l")
    plt.figure()
    plt.plot(af.drag_coefficient, af.lift_coefficient)
    plt.grid()
    plt.xlabel("C_d"), plt.ylabel("C_l")
    plt.show() 
    pass 
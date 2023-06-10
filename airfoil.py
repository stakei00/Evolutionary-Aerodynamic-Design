import numpy as np 
import math 

class Airfoil:

    def __init__(self, m:float, p:float, t:float, Re:int, aseq:list=[-15,24,1],\
                 n_iter:int=300, ncrit:int=9, airfoil_hist:dict=None) -> None:
        """
        Initializes airfoil object
        Inputs:
            m = max camber as percentage of chord 
            p = max camber location 
            t = thickness to chord ratio 
            Re = Reynold's number of airfoil (chord = 1) in increments of 100,000 
            airfoil_stats: dictionary containing solutions to previously run airfoils 
        """
        self.max_camber = m
        self.max_camber_loc = p
        self.thickness_to_chord = t
        self.Re = Re
        
        #create equivalent NACA 4-digit designation
        if round(m*100, 0) < 1: 
            self.NACA_4series_desig = "00"
        else: 
            self.NACA_4series_desig = str(int(10*round(m*100, 0) + round(p*10, 0)))
        if round(t*100,0) < 10: 
            self.NACA_4series_desig += "0"
        self.NACA_4series_desig += str(int(round(t*100, 0)))                     

        #if airfoil history dictionary is provided, see if self is in it 
        if airfoil_hist is not None: 
            if self.NACA_4series_desig in airfoil_hist.keys():
                self.CDCL_avl = airfoil_hist[self.NACA_4series_desig][0]
                self.max_lift_coefficient = airfoil_hist[self.NACA_4series_desig][1]
                self.Re = airfoil_hist[self.NACA_4series_desig][2]
                return 

        #if self not found in the history, run XFOIL to get attributes
        xfoil_data_neg = self.run_XFOIL_aseq(0, aseq[0], aseq[2], self.Re,\
                            n_iter=n_iter, ncrit=ncrit) #negative angles of attack 
        xfoil_data_pos = self.run_XFOIL_aseq(1, aseq[1], aseq[2], self.Re,\
                            n_iter=n_iter, ncrit=ncrit) #positive angles of attack
        
        if None in [xfoil_data_neg, xfoil_data_pos]:
            return

        xfoil_data_neg = np.flipud(xfoil_data_neg)
        xfoil_data = np.concatenate((xfoil_data_neg, xfoil_data_pos))

        self.alpha = xfoil_data[:,0]
        self.lift_coefficient = xfoil_data[:,1]
        self.drag_coefficient = xfoil_data[:,2]
        self.moment_coefficient = xfoil_data[:,4]

        self.find_3pt_drag_polar()
        
    def run_XFOIL_aseq(self, alpha_i:float, alpha_f:float, alpha_step:float, \
                  Re:int or float, n_iter:int, ncrit:int) -> None: 
        """
        Analysis of airfoil in XFOIL. Adds results to object attributes
        """
        import os
        import subprocess
        import os
        import time

        if os.path.exists("polar_file.txt"):
            os.remove("polar_file.txt")

        input_content = "PLOP\n" +\
            "G F\n" +\
            "\n" +\
            f"NACA {self.NACA_4series_desig}\n" +\
            "PANE\n" +\
            "OPER\n" +\
            "VPAR\n" +\
            "N\n" +\
            f"{ncrit}\n\n" +\
            f"Visc {Re}\n" +\
            "PACC\n" +\
            "polar_file.txt\n\n" +\
            f"ITER {n_iter}\n" +\
            f"ASeq {alpha_i} {alpha_f} {alpha_step}\n" +\
            "\n\n" +\
            "quit\n" 
                
        timeout_seconds = 10 #process will terminate after this amount of time 

        with open(os.devnull, 'w') as null_file: 
            process = subprocess.Popen(["xfoil.exe"], stdin=subprocess.PIPE, stdout=null_file)      
            process.stdin.write(input_content.encode())
            process.stdin.close()

        start_time = time.time()
        while time.time() - start_time < timeout_seconds:
            if process.poll() is not None: #check if the process has completed
                break
            time.sleep(0.1) 

        if process.poll() is None:
            process.terminate()
            return None 

        try: 
            return np.loadtxt("polar_file.txt", skiprows=12)
        except: 
            return None
        
    def find_3pt_drag_polar(self) -> None: 
        """
        Gets 3 points on CLCD drag polar for use with AVL
        """
        alpha_interp = np.linspace(min(self.alpha), max(self.alpha), 200)
        cl_interp = np.interp(alpha_interp, self.alpha, self.lift_coefficient)
        cl_max_ind = np.argmax(cl_interp)
        cl_min_ind = np.argmin(cl_interp)
        
        #trim alpha, cl, cd arrays 
        alpha_new, cl_new, cd_new, cm_new = [],[],[],[]
        for i,a in enumerate(self.alpha):
            if a < alpha_interp[cl_min_ind]: continue
            if a > alpha_interp[cl_max_ind]: continue

            alpha_new.append(a)
            cm_new.append(self.moment_coefficient[i])
            cl_new.append(self.lift_coefficient[i])
            cd_new.append(self.drag_coefficient[i]) 

        #update attribute arrays 
        self.alpha = alpha_new 
        self.lift_coefficient = cl_new
        self.drag_coefficient = cd_new
        self.moment_coefficient = cm_new 

        #interpolate drag coefficient
        if len(self.lift_coefficient) < 2: 
            return 
        cd_interp = np.interp(cl_interp, self.lift_coefficient, self.drag_coefficient)
        cd_min_ind = np.argmin(cd_interp)

        #pack up everything and add as object attribute
        cl1, cd1 = cl_interp[cl_min_ind], cd_interp[cl_min_ind]
        cl2, cd2 = cl_interp[cd_min_ind], cd_interp[cd_min_ind]
        cl3, cd3 = cl_interp[cl_max_ind], cd_interp[cl_max_ind]

        self.CDCL_avl = [cl1, cd1, cl2, cd2, cl3, cd3]
        self.max_lift_coefficient = cl_interp[cl_max_ind] 

    def try_add_to_hist(self, airfoil_hist:dict) -> dict:
        
        if airfoil_hist is None or not hasattr(self, "CDCL_avl"): 
            return 
        
        if self.NACA_4series_desig not in airfoil_hist.keys(): 
            airfoil_hist[self.NACA_4series_desig] = [self.CDCL_avl, \
                                                  self.max_lift_coefficient, self.Re]
        return airfoil_hist

    def plot_airfoil(self, ax, chord_scale:float=None, origin:tuple or list=None,\
                      twist_deg:float=None, linecolor:str="black") -> None:
        """
        plots naca 4 series airfoil
        """
        m,p,t = self.max_camber, self.max_camber_loc, self.thickness_to_chord
        x_array = np.linspace(0,1,200)
        yt_array = []
        yc_array = []
        thet_array = []
        for x in x_array:
            yt_array.append(5*t*(0.2929*math.sqrt(x) - 0.126*x - 0.3516*x**2 + \
                                 0.2843*x**3 - 0.1015*x**4))
            if 0<=x<=p: 
                yc_array.append(m*(2*p*x - x**2)/(p**2))
                dycdx = 2*m*(p-x)/(p**2)
            elif p<x<=1:
                yc_array.append(m*((1-2*p) + 2*p*x - x**2)/((1-p)**2))
                dycdx = 2*m*(p-x)/((1-p)**2)

            thet_array.append(math.atan(dycdx))

        x_upper = np.array([x - yt_array[i]*math.sin(thet_array[i]) \
                            for i,x in enumerate(x_array)])
        x_lower = np.array([x + yt_array[i]*math.sin(thet_array[i]) \
                            for i,x in enumerate(x_array)])
        x = np.concatenate((x_upper, np.flip(x_lower)))
        y_upper = np.array([yc + yt_array[i]*math.cos(thet_array[i]) \
                            for i,yc in enumerate(yc_array)])
        y_lower = np.array([yc - yt_array[i]*math.cos(thet_array[i]) \
                            for i,yc in enumerate(yc_array)])
        y = np.concatenate((y_upper, np.flip(y_lower))) 

        #scale airfoil based on chord
        if chord_scale is not None: 
            x = np.multiply(-chord_scale, x)
            y = np.multiply(chord_scale, y)

        #twist 
        if twist_deg is not None: 
            twist_rad = math.radians(-twist_deg)
            rotation_matrix = np.array([[np.cos(twist_rad), np.sin(twist_rad)],
                                    [-np.sin(twist_rad), np.cos(twist_rad)]])
            xys = np.vstack((x,y))
            x,y = np.dot(rotation_matrix, xys)

        #add offset
        if origin is not None:
            x = np.add(origin[0], x)
            y = np.add(origin[1], y)

        #plot them lines! 
        ax.plot(x, y, color=linecolor, linewidth=1)
        

if __name__ == "__main__":
    m = 0.02
    p = 0.4
    t = 0.12
    re = 10e6

    af = Airfoil(m,p,t,re)

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

import numpy as np 
import math 

class Airfoil:

    def __init__(self, m:float, p:float, t:float) -> None:
        """
        Initializes airfoil object 
        """
        self.max_camber = m
        self.max_camber_loc = p
        self.thickness_to_chord = t
        
        self.NACA_4series_desig = str(int(1000*round(m*100, 0) + 100*round(p*10, 0) + round(t*100, 0)))
        self.run_XFOIL(-12, 15, 1, 10e6, 50) #test numbers 
        self.find_3pt_drag_polar()

    def run_XFOIL(self, alpha_i:float, alpha_f:float, alpha_step:float, \
                  Re:int or float, n_iter:int) -> None: 
        """
        Analysis of airfoil in XFOIL. Adds results to object attributes
        """
        import os
        import subprocess
        import numpy as np

        if os.path.exists("polar_file.txt"):
            os.remove("polar_file.txt")

        input_file = open("input_file.in", 'w')
        input_file.write("PLOP\n")
        input_file.write("G F\n")
        input_file.write("\n")
        input_file.write(f"NACA {self.NACA_4series_desig}\n")
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

    def find_3pt_drag_polar(self) -> None: 
        """
        Gets 3 points on CLCD drag polar for use with AVL 
        """
        #get minimum drag coefficient and corresponding lift coefficient: 
        i_Cdmin = np.argmin(self.drag_coefficient)
        Cl2, Cd2 = self.lift_coefficient[i_Cdmin], self.drag_coefficient[i_Cdmin]

        #find stall Cl&Cd for positive and negative Cl regions 
        h = math.radians(self.alpha[1] - self.alpha[0]) 
        first_derivative = np.gradient(self.lift_coefficient,h)
        second_derivative = np.gradient(first_derivative, h) 

        Cl1, Cd1 = self.lift_coefficient[0], self.drag_coefficient[0]
        

        #negative Cl region
        for i,Cl in enumerate(self.lift_coefficient): 
            if Cl > 0: continue
            if abs(second_derivative[i]) > 10: 
                Cl1, Cd1 = self.lift_coefficient[i], self.drag_coefficient[i]

        #positive Cl region
        for i,Cl in enumerate(self.lift_coefficient):
            if Cl < 0: continue
            if abs(second_derivative[i]) > 10:  
                Cl3, Cd3 = self.lift_coefficient[i], self.drag_coefficient[i]

        self.CDCL_avl = [Cl1, Cd1, Cl2, Cd2, Cl3, Cd3]
        self.max_lift_coefficient = Cl3

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
            yt_array.append(5*t*(0.2929*math.sqrt(x) - 0.126*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4))
            if 0<=x<=p: 
                yc_array.append(m*(2*p*x - x**2)/(p**2))
                dycdx = 2*m*(p-x)/(p**2)
            elif p<x<=1:
                yc_array.append(m*((1-2*p) + 2*p*x - x**2)/((1-p)**2))
                dycdx = 2*m*(p-x)/((1-p)**2)

            thet_array.append(math.atan(dycdx))

        x_upper = np.array([x - yt_array[i]*math.sin(thet_array[i]) for i,x in enumerate(x_array)])
        x_lower = np.array([x + yt_array[i]*math.sin(thet_array[i]) for i,x in enumerate(x_array)])
        x = np.concatenate((x_upper, np.flip(x_lower)))
        y_upper = np.array([yc + yt_array[i]*math.cos(thet_array[i]) for i,yc in enumerate(yc_array)])
        y_lower = np.array([yc - yt_array[i]*math.cos(thet_array[i]) for i,yc in enumerate(yc_array)])
        y = np.concatenate((y_upper, np.flip(y_lower))) 

        #scale airfoil based on chord
        if chord_scale is not None: 
            x = np.multiply(-chord_scale, x)
            y = np.multiply(chord_scale, y)

        #twist 
        if twist_deg is not None: 
            twist_rad = math.radians(twist_deg)
            rotation_matrix = np.array([[np.cos(twist_rad), np.sin(twist_rad)],
                                    [-np.sin(twist_rad), np.cos(twist_rad)]])
            xys = np.vstack((x,y))
            x,y = np.dot(rotation_matrix, xys)

        #add offset
        if origin is not None:
            x = np.add(origin[0], x)
            y = np.add(origin[1], y)

        #plot them lines! 
        ax.plot(x, y, color=linecolor)
        

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
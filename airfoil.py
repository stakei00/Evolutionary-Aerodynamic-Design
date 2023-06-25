import numpy as np
from scipy.interpolate import interp1d 
import math 
import pickle
import os

class Airfoil:

    def __init__(self, m:float, p:float, t:float, Re:int, aseq:list=[-24,24,1],\
                 n_iter:int=200, ncrit:int=9, search_airfoils=False) -> None:
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
        self.Re = int(round(Re, -5))
        
        #create equivalent NACA 4-digit designation
        if round(m*100, 0) < 1 or round(p*10,0) < 1: 
            self.NACA_4series_desig = "00"
        else: 
            self.NACA_4series_desig = str(int(10*round(m*100, 0) + round(p*10, 0)))
        if round(t*100,0) < 10: 
            self.NACA_4series_desig += "0"
        self.NACA_4series_desig += str(int(round(t*100, 0)))

        #check if airfoil already exists in pickle depot
        filepath = f"airfoils/n{self.NACA_4series_desig}_{self.Re}"
        if os.path.exists(filepath) and search_airfoils:
            with open(filepath, "rb") as file:
                airfoil = pickle.load(file)
                self.alpha = airfoil.alpha
                self.lift_coefficient = airfoil.lift_coefficient
                self.drag_coefficient = airfoil.drag_coefficient 
                self.moment_coefficient = airfoil.moment_coefficient
                self.CDCL_avl = airfoil.CDCL_avl
                self.max_lift_coefficient = airfoil.max_lift_coefficient
            return 

        #if self not found in the history, run XFOIL to get attributes
        xfoil_data_neg = self.run_XFOIL_aseq(0, aseq[0], aseq[2], self.Re,\
                            n_iter=n_iter, ncrit=ncrit) #negative angles of attack 
        xfoil_data_pos = self.run_XFOIL_aseq(1, aseq[1], aseq[2], self.Re,\
                            n_iter=n_iter, ncrit=ncrit) #positive angles of attack

        if np.size(xfoil_data_neg) == 7: 
            xfoil_data_neg = np.reshape(xfoil_data_neg, (1,7))
        if np.size(xfoil_data_pos) == 7:
            xfoil_data_pos = np.reshape(xfoil_data_pos, (1,7))
        try: 
            if np.size(xfoil_data_neg) > 7:
                xfoil_data_neg = np.flipud(xfoil_data_neg)
            xfoil_data = np.concatenate((xfoil_data_neg, xfoil_data_pos))
            self.alpha_raw = xfoil_data[:,0]
            self.lift_coefficient_raw = xfoil_data[:,1]
            self.drag_coefficient_raw = xfoil_data[:,2]
            self.moment_coefficient_raw = xfoil_data[:,4]
        except: 
            return

        self.find_3pt_drag_polar()
        
        if search_airfoils:
            self.pickle_airfoil()
        
    def run_XFOIL_aseq(self, alpha_i:float, alpha_f:float, alpha_step:float, \
                  Re:int or float, n_iter:int, ncrit:int) -> None: 
        """
        Analysis of airfoil in XFOIL. Adds results to object attributes
        """
        import os
        import subprocess
        import os
        import time
        import random 

        #generate unique file name for polar data 
        while True: 
            polar_tag = str(int(random.uniform(0, 100000)))
            polar_file = f"polar_{polar_tag}.txt"
            if not os.path.exists(polar_file):
                break 

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
            f"{polar_file}\n\n" +\
            f"ITER {n_iter}\n" +\
            f"ASeq {alpha_i} {alpha_f} {alpha_step}\n" +\
            "\n\n" +\
            "quit\n" 
                
        timeout_seconds = 12 #process will terminate after this amount of time 

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
            print(f"\nXFOIL process on NACA{self.NACA_4series_desig} timed out." +\
                  f" Re = {self.Re}. Returning incomplete polar")
            process.terminate()

        try: 
            polar = np.loadtxt(polar_file, skiprows=12)
            os.remove(polar_file)
            return polar 
        
        except: 
            return None
        
    def find_3pt_drag_polar(self) -> None: 
        """
        Gets 3 points on CLCD drag polar for use with AVL
        """

        #expand raw data points using interpolation
        self.alpha_interp = np.linspace(min(self.alpha_raw), max(self.alpha_raw), 200)
        cl_interp_func = interp1d(self.alpha_raw, self.lift_coefficient_raw)
        cd_interp_func = interp1d(self.alpha_raw, self.drag_coefficient_raw)
        self.cl_interp = np.array([cl_interp_func(a) for a in self.alpha_interp])
        dcl_dalpha = np.gradient(self.cl_interp, self.alpha_interp[1]-self.alpha_interp[0]) #lift slope

        #get index of first angle of attack >= 0: 
        ind_alpha0 = int(len(self.alpha_interp)/2)
        for i,a in enumerate(self.alpha_interp[1:]):
            if a >= 0 and self.alpha_interp[i-1] < 0:
                ind_alpha0 = i
                break 

        baseline_lift_slope = dcl_dalpha[ind_alpha0]

        #find max lift coefficient
        cl_max_ind = None
        for i,a in enumerate(self.alpha_interp[:-1]):
            if a < 0: continue 
            if dcl_dalpha[i]*dcl_dalpha[i+1] < 0: #if sign change 
                cl_max_ind = i
                break 
                
        if cl_max_ind is None: 
            cl3, cd3 = None, None
            self.max_lift_coefficient = None
        else:     
            self.alpha_max_lift_coeff = self.alpha_interp[cl_max_ind]
            self.max_lift_coefficient = self.cl_interp[cl_max_ind]
            cl3, cd3 = self.max_lift_coefficient, cd_interp_func(self.alpha_max_lift_coeff)

        #try to find min lift coefficient
        alpha_interp_flip = np.flip(self.alpha_interp)
        dcl_dalpha_flip = np.flip(dcl_dalpha)

        cl_min_ind = None
        for i,a in enumerate(alpha_interp_flip[:-1]):
            if a > 0: continue 
            if dcl_dalpha_flip[i]*dcl_dalpha_flip[i+1] < 0: #if sign change 
                cl_min_ind = len(self.alpha_interp)-1-i 
                break

        if cl_min_ind is None: 
            cl_min_ind = 0
            cl1, cd1 = None, None
            self.min_lift_coefficient = None 
        else:  
            self.min_lift_coefficient = self.cl_interp[cl_min_ind]
            self.alpha_min_lift_coeff =  self.alpha_interp[cl_min_ind]
            cl1, cd1 = self.min_lift_coefficient, cd_interp_func(self.alpha_min_lift_coeff)
                
        #trim alpha, cl, cd arrays 
        alpha_new, cl_new, cd_new, cm_new = [],[],[],[]
        for i,a in enumerate(self.alpha_raw):
            if cl_min_ind is not None: 
                if a < self.alpha_interp[cl_min_ind]: continue
            if cl_max_ind is not None: 
                if a > self.alpha_interp[cl_max_ind]: break

            alpha_new.append(a)
            cm_new.append(self.moment_coefficient_raw[i])
            cl_new.append(self.lift_coefficient_raw[i])
            cd_new.append(self.drag_coefficient_raw[i])
        
        self.alpha = alpha_new 
        self.lift_coefficient = cl_new
        self.moment_coefficient = cm_new
        self.drag_coefficient = cd_new

        #interpolate drag coefficient
        if len(self.lift_coefficient) < 2: 
            return
        
        self.alpha_interp = self.alpha_interp[cl_min_ind:cl_max_ind+1]
        self.cl_interp = np.array([cl_interp_func(a) for a in self.alpha_interp])
        self.cd_interp = np.array([cd_interp_func(a) for a in self.alpha_interp])
        cd_min_ind = np.argmin(self.cd_interp)

        dcl_dalpha = dcl_dalpha[cl_min_ind:cl_max_ind+1]
        for i,a in enumerate(self.alpha_interp):
            if a < 0: continue
            if dcl_dalpha[i] < 0.2*baseline_lift_slope:
                cl3 = self.cl_interp[i]
                cd3 = self.cd_interp[i]
                break 
        
        alpha_interp_flip = np.flip(self.alpha_interp)
        dcl_dalpha_flip = np.flip(dcl_dalpha)
        for i,a in enumerate(alpha_interp_flip):
            if a >= 0: continue 
            if dcl_dalpha_flip[i] < 0.2*baseline_lift_slope: 
                ind = len(self.cl_interp)-1-i
                cl1 = self.cl_interp[ind]
                cd1 = self.cd_interp[ind]
        
        #pack up everything and add as object attribute
        cl2, cd2 = self.cl_interp[cd_min_ind], self.cd_interp[cd_min_ind]

        self.CDCL_avl = [cl1, cd1, cl2, cd2, cl3, cd3]

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
        line, = ax.plot(x, y, color=linecolor, linewidth=1)
        return line 

    def pickle_airfoil(self):
        """
        saves the airfoil to pickle
        """

        with open(f"airfoils/n{self.NACA_4series_desig}_{self.Re}", "wb") as file:
            pickle.dump(self, file)

class Airfoil_Interpolator:
    """
    airfoil data interpreter. Currently designed for naca 4 digit airfoils with 
    independent parameters: (m,p,t,Re)
    """ 

    def __init__(self, batch_data_file:str):
        
        self.batch_data_file_path = batch_data_file
        #load batch data 
        self.load_batch_data()
        #create interpolator from batch data
        self.create_interpolator()

    def load_batch_data(self) -> None:
        
        import json
        self.batch_data = json.load(open(self.batch_data_file_path, "r"))
        
        #create point range
        m_list, p_list, t_list, re_list = [],[],[],[]
        for case in self.batch_data["cases"]: 
            #get m,p,t from airfoil naca digits
            m_list.append(case[0])
            p_list.append(case[1])
            t_list.append(case[2])
            re_list.append(case[3])

        #get rid of duplicate values 
        m_list = [*set(m_list)]
        p_list = [*set(p_list)]
        t_list = [*set(t_list)]
        re_list = [*set(re_list)]

        #sort lists to get ranges
        m_list.sort()
        p_list.sort() 
        t_list.sort() 
        re_list.sort() 

        points = (m_list, p_list, t_list, re_list)
        values = np.zeros((len(m_list), len(p_list), len(t_list), len(re_list), 7)) 

        for i,m in enumerate(m_list):
            for j,p in enumerate(p_list): 
                for k,t in enumerate(t_list): 
                    for l,re in enumerate(re_list): 
                        
                        naca4 = str(int(10*round(m*100, 0) + round(p*10, 0)))
                        if round(m*100, 0) < 1: 
                            naca4 = "00"
                        if round(t*100,0) < 10: 
                            naca4 += "0"
                        naca4 += str(int(round(t*100, 0)))

                        key = naca4+"_"+str(int(re))

                        if key not in list(self.batch_data.keys()):
                            values[i,j,k,l] = 7*[math.nan]
                        else: 
                            data = self.batch_data[key]
                            for ii,val in enumerate(data):
                                if val == "None": 
                                    data[ii] = math.nan 
                                    break

                            values[i,j,k,l] = [data[0],data[2][0],\
                                               data[2][1],data[2][2],\
                                                data[2][3],data[2][4],data[2][5]]
                         
        self.points = points 
        self.values = values 

    def create_interpolator(self):

        import scipy.interpolate as sci_int

        self.interpolator = sci_int.RegularGridInterpolator(self.points,\
                                self.values, bounds_error=False,\
                                    fill_value=None)


if __name__ == "__main__":
    
    #testing out custom airfoil 
    afoil = Airfoil(m=0.06, p=0.3, t=0.09, Re=1e6)


    #testing out airfoil interpolator
    interpolator = Airfoil_Interpolator("xfoil_batch_data.json")

    test = interpolator.interpolator([0.04, 0.2, 0.12, 300000])
    print("airfoil NACA4412 at Re=300000")
    print(f"CDCL input: {test[0][1:]}")
    print(f"CL max: {test[0][0]}")

    pass 
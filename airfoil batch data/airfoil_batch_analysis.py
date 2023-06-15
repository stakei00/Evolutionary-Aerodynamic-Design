import sys, os 
sys.path.append(os.getcwd())
import airfoil
import matplotlib.pyplot as plt
import json
import random

"""
runs XFOIL for a massive batch of NACA 4 series airfoils. Dumps to json file. 
Preliminary step for interpolator creation. 
"""

def analyze_airfoil(inputs):
    [m,p,t,re] = inputs
    try:
        af = airfoil.Airfoil(m, p, t, re, aseq=[-24,24,1],\
                 n_iter=500, ncrit=9, search_airfoils=False)
     
        res = {}
        min_lift_coefficient = "None"

        if hasattr(af, min_lift_coefficient):
            min_lift_coefficient = af.min_lift_coefficient

        key = af.NACA_4series_desig+"_"+str(af.Re)
        res[key] = [af.max_lift_coefficient, min_lift_coefficient, af.CDCL_avl]
      
        return res
    except: return None 

if __name__ == "__main__": 

    m_list = [x/100 for x in list(range(0, 9))]
    p_list = [x/10 for x in list(range(0, 9))]
    t_list = [x/100 for x in list(range(7, 25))]
    re_list = [x*1e5 for x in list(range(1,11))] 

    #get permutations
    try: 
        batch_res_dict = json.load(open("xfoil_batch_data.json", "r"))
        tried = batch_res_dict["cases"]
    except: 
        batch_res_dict = {}
        tried = []
    cases = []
    timeout = []


    for m in m_list:
        for p in p_list:
            for t in t_list: 
                for re in re_list: 
                    cases.append([m,p,t,re])

    plt.ion()
    fig, axs = plt.subplots(2,2,figsize=(10,10))


    def calculate():
        
        while True: 
            case = random.choice(cases)
            if case not in tried and case not in timeout:
                break 
        print
        res = analyze_airfoil(case)
        if res is None: 
            timeout.append(case)
            return None
        
        tried.append(case)
        
        #update the plot
        axs[0,0].clear(), axs[0,1].clear(), axs[1,0].clear(), axs[1,1].clear()
        axs[0,0].scatter([lis[0] for lis in tried], [lis[1] for lis in tried], c=[lis[3] for lis in tried])
        axs[0,1].scatter([lis[0] for lis in tried], [lis[2] for lis in tried], c=[lis[3] for lis in tried])
        axs[1,0].scatter([lis[1] for lis in tried], [lis[2] for lis in tried], c=[lis[3] for lis in tried])
        axs[1,1].scatter([lis[0] for lis in tried], [lis[1] for lis in tried], c=[lis[2] for lis in tried])

       
        fig.canvas.draw()
        fig.canvas.flush_events()
        
        return res
    
    af_dicts = []
    count = 1
    tot_count = 500 #! number of iterations to run through
    while count < tot_count: 
        
        print(f"\tsolving airfoil: {count}/{tot_count}")
        af = calculate()
        if af is not None: 
            af_dicts.append(af)
        count += 1
        #pool = mp.Pool(processes=mp.cpu_count())
        #af_dicts = pool.map(calculate, [])

    with open("xfoil_batch_data.json", "w") as file: 
        for dic in af_dicts:  
            batch_res_dict[list(dic.keys())[0]] = dic[list(dic.keys())[0]]

        batch_res_dict["cases"] = tried
        json.dump(batch_res_dict, file)

    pass 
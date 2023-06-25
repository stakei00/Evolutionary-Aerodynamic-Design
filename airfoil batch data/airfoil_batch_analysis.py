import sys, os 
sys.path.append(os.getcwd())
import airfoil
import matplotlib.pyplot as plt
import pandas as pd
import random
import multiprocessing as mp

"""
runs XFOIL for a massive batch of NACA 4 series airfoils. Dumps to json file. 
Preliminary step for interpolator creation. 
"""

def analyze_airfoil(inputs):
    [m,p,t,re] = inputs
    
    af = airfoil.Airfoil(m, p, t, re, aseq=[-24,24,1],\
             n_iter=500, ncrit=9, search_airfoils=False)
    try: 
        res = af.CDCL_avl + [af.max_lift_coefficient, af.min_lift_coefficient]
        print(f"\t{[inputs] + [round(float(x),3) for x in res if x is not None]}")
    except: 
        return [None]*8
    if res is None: 
        return [None]*8
    return res


def calculate(inputs):
        cases, tried = inputs
        while True: 
            case = random.choice(cases)
            if case[0] == 0 or case[1] == 0:
                case[:2] = [0,0]

            if case not in tried:
                break 

        af = analyze_airfoil(case)

        dict_ = {
            "m":        case[0],
            "p":        case[1],
            "t":        case[2],
            "Re":       case[3],
            "cl1":      af[0],
            "cd1":      af[1],
            "cl2":      af[2],
            "cd2":      af[3],
            "cl3":      af[4],
            "cd3":      af[5],
            "clmax":    af[6],
            "clmin":    af[7]
            }
        
        af_df = pd.DataFrame(data=dict_, index=[0])
        return [af_df,case]


if __name__ == "__main__": 

    m_list = [x/100 for x in list(range(0, 9))]
    p_list = [x/10 for x in list(range(0, 9))]
    t_list = [x/100 for x in list(range(7, 25))]
    re_list = [x*1e5 for x in list(range(1,11))] 

    batch_res_df = pd.read_csv("airfoil batch data/xfoil_batch_data.csv")
    initial_count = batch_res_df.shape[0]
    cases = []
    for m in m_list:
        for p in p_list:
            for t in t_list: 
                for re in re_list: 
                    cases.append([m,p,t,re])

    tried = [row[:4] for row in batch_res_df.values.tolist()]
    tot_count = 100 #! number of iterations to run through
    num_workers = mp.cpu_count()-2
    with mp.Pool(num_workers) as pool: 
        while batch_res_df.shape[0]-initial_count <= tot_count: 
            
            print(f"\n\t\tairfoils: {batch_res_df.shape[0]-initial_count}/{tot_count}\n")
            try: results = pool.map(calculate, [(cases,tried)]*num_workers)
            except: continue
            for res in results:
                df, case = res 
                batch_res_df = pd.concat([batch_res_df, df], ignore_index=True)
                tried.append(case)

    batch_res_df = batch_res_df.sort_values(['m','p','t','Re'])
    batch_res_df.to_csv("airfoil batch data/xfoil_batch_data.csv", index=False)
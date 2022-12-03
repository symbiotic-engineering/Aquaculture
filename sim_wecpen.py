import modules
from modules import Aqua_Obj
from optimization import OpData
import optimization
from utilities import *
import copy
import numpy as np
import importlib
importlib.reload(modules)
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

max_iter = 2000

def wecpen_opt(all_vars, *x0):
    wec_types = ['point_absorber_RM3'] #,'attenuator','terminator','point_absorber_RM3']

    # design variables and parameters
    if 'p_wave_vec' in all_vars:
        x_name = ['x_wec','x_pen', 'x_es']
        p_name = ['x_type_wec', 'x_env'] #,'x_env'
    else: # es is not required when average wave is given
        x_name = ['x_wec','x_pen']
        p_name = ['x_type_wec', 'x_env', 'x_es'] #,'x_env'

    x = OpData(x_name)
    if x0:
        for i in range(len(x.list)):
            x.nom_dict[x.list[i]] = x0[0][i]

    x_disc_name = ['x_disc_pen']
    x_disc = OpData(x_disc_name)
    #x_disc.bnds[0] = (11, 11)   # to accelerate the computation

    
    
    param = OpData(p_name + x_disc_name)

    #Environemental parameter can be set either here or in module.
    #param.nom_dict['temp'] =     16  #'C'
    #param.nom_dict['salinity'] = 33  #'[PSU]'
    #param.nom_dict['U'] =        0.2 #'[m/s]'
    #param.nom_dict['O2_in'] =    8   #'[mg/l]'
    #param.nom_dict['wave_height'] = 1.4  #'[m]'     
    #param.nom_dict['wave_period'] = 8.33 #'[s]'
    #param_val['wave_data'] = "wave_data/Wave_Data2.csv"

    #optimization
    res={}
    res_best = {}
    init_flag = 1

    for i in range(len(wec_types)):
        param.nom_dict['wec_type'] = wec_types[i]
        for j0 in range(len(x_disc.list)):
            lower_bnd, upper_bnd = x_disc.bnds[j0]
            for j1 in range(lower_bnd,upper_bnd+1):
                param.nom_dict[x_disc.list[j0]] = j1            

                if 'x_wave_ave' in all_vars:
                    param.nom_dict['es_size'] = 0                    

                res_opt, op_obj, p = optimization.run_optimization(x.name, x.nom0, param.name, param.nom_dict, all_vars, max_iter)
                if init_flag:
                    x_init, p_init = x , p
                    res_best = copy.copy(res_opt)
                    p_best = p.nom_dict
                elif (res_opt.success) and (res_opt.fun < res_best.fun):
                    res_best = copy.copy(res_opt)
                    p_best = p.nom_dict

                init_flag = 0
    
    return x_init, p_init, x, res_best, p_best, op_obj

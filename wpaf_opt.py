from optimization import OpData
import optimization
from utilities import *
import copy
import numpy as np


max_iter = 2000

def wpaf_opt(all_vars, *args):
    wec_types = ['point absorber'] #,'attenuator','terminator','point_absorber_RM3']

    # design variables and parameters
    if 'p_wave_vec' in all_vars:
        x_name = ['x_wec','x_pen', 'x_es']
        p_name = ['x_type_wec', 'x_env', 'gis_handler'] #,'x_env'
    else: # es is not required when average wave is given
        x_name = ['x_wec','x_pen']
        p_name = ['x_type_wec', 'x_env', 'x_es', 'gis_handler'] #,'x_env'

    x = OpData(x_name)
    if 'x0' in args[0]:
        for i in range(len(x.list)):
            x.nom_dict[x.list[i]] = args[0]['x0'][i]

    x_disc_name = ['x_disc_pen']
    x_disc = OpData(x_disc_name)
    if 'fixed_num_pen' in args[0]:
        x_disc.bnds[0] = (args[0]['fixed_num_pen'], args[0]['fixed_num_pen'])   # to accelerate the computation
    
    
    param = OpData(p_name + x_disc_name)

    if 'handler' in args[0]:
        param.nom_dict['handler'] = args[0]['handler']
    else:
        print('GIS handler is needed!')
        exit()
    
    #optimization
    res = {}
    res_best = {}
    init_flag = 1

    for i in range(len(wec_types)):
        param.nom_dict['wec_type'] = wec_types[i]

        for x_disc_i in range(len(x_disc.list)):
            lower_bnd, upper_bnd = x_disc.bnds[x_disc_i]

            for x_disc_val in range(lower_bnd,upper_bnd+1):
                param.nom_dict[x_disc.list[x_disc_i]] = x_disc_val            

                if 'x_wave_ave' in all_vars:
                    param.nom_dict['es_size'] = 0                    

                print('optimization is running ...')
                res, op_obj, p = optimization.run_optimization(x.name, x.nom0, param.name, param.nom_dict, all_vars, max_iter)

                if init_flag:
                    x_init, p_init = x, p
                    res_best = copy.deepcopy(res)
                    p_best = p.nom_dict
                elif (res.success) and (res.fun < res_best.fun):
                    res_best = copy.deepcopy(res)
                    p_best = p.nom_dict

                init_flag = 0
    
    return x_init, p_init, x, res_best, p_best, op_obj

# ============================================================================ #
#                       Multi-Objective Optimization                           #
# ============================================================================ #

def wpaf_multi_opt(all_vars, *args):
    wec_types = ['point absorber']
    
    # design variables and parameters
    if 'p_wave_vec' in all_vars:
        x_name = ['x_wec','x_pen', 'x_es']
        p_name = ['x_type_wec', 'x_env', 'gis_handler'] #,'x_env'
    else: # es is not required when average wave is given
        x_name = ['x_wec','x_pen']
        p_name = ['x_type_wec', 'x_env', 'x_es', 'gis_handler'] #,'x_env'
        
    x = OpData(x_name)
    if 'x0' in args[0]:
        for i in range(len(x.list)):
            x.nom_dict[x.list[i]] = args[0]['x0'][i]

    x_disc_name = ['x_disc_pen']
    x_disc = OpData(x_disc_name)
    if 'fixed_num_pen' in args[0]:
        x_disc.bnds[0] = (args[0]['fixed_num_pen'], args[0]['fixed_num_pen'])
    
    param = OpData(p_name + x_disc_name)

    if 'handler' in args[0]:
        param.nom_dict['handler'] = args[0]['handler']
    else:
        print('GIS handler is needed!')
        exit()
    
    result = []
    
    for i in range(len(wec_types)):
        param.nom_dict['wec_type'] = wec_types[i]

        for x_disc_i in range(len(x_disc.list)):
            lower_bnd, upper_bnd = x_disc.bnds[x_disc_i]

            for x_disc_val in range(lower_bnd,upper_bnd+1):
                param.nom_dict[x_disc.list[x_disc_i]] = x_disc_val            

                if 'x_wave_ave' in all_vars:
                    param.nom_dict['es_size'] = 0                    
                
                print('multi optimization is running ...')
    
                des_var, op_obj = optimization.run_multi_optimization(x.name, param.name, param.nom_dict, all_vars, max_iter)
                
    result.append((des_var, op_obj))
    
    return result

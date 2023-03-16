from optimization import OpData
import optimization
from utilities import *
import copy
import numpy as np

max_iter = 2000    

def env_opt(all_vars, *args):
    wec_types = ['point absorber'] #,'attenuator','terminator','point_absorber_RM3']

    # design variables
    x_name = ['pos_env']
    x = OpData(x_name)

    # parameters
    param_name = ['x_wec','x_type_wec', 'x_pen', 'gis_handler']
    param = OpData(param_name)

    # WEC and Pen Parameters are defined by optimal results obtained by running "run_sim_wec" 
    param.nom_dict['capture_width']=  np.NaN   #[m]
    param.nom_dict['pen_diameter']=   21.4     #[m]
    param.nom_dict['pen_height']=     19.3     #[m]
    param.nom_dict['stock_density']=  20       #[kg/m^3]

    if 'handler' in args[0]:
        param.nom_dict['handler'] = args[0]['handler']
    else:
        print('GIS handler is needed!')
        exit()

    #optimization
    res={}
    res_best={}
    init_flag = 1

    for i in range(len(wec_types)):
        param.nom_dict['wec_type'] = wec_types[i]
        
        res, op_obj, p = optimization.run_optimization(x.name, x.nom0, param.name, param.nom_dict, all_vars, max_iter)
        
        if init_flag:
            x_init, p_init = x , p
            res_best = copy.deepcopy(res)
            p_best = p.nom_dict
        elif (res.success) and (res.fun < res_best.fun):
            res_best = copy.deepcopy(res)
            p_best = p.nom_dict
        
        init_flag = 0
    
    return x_init, p_init, x, res_best, p_best, op_obj  
from gis.gis_handler import GISHandler
from optimization import OpData
import optimization
from utilities import *
import copy
import numpy as np
import random
#import importlib
#importlib.reload(modules)

max_iter = 2000

all_vars_default = ['x_wec','x_type_wec','x_pen','p_pen','x_env','p_env','p_wec','p_fish_salmon','pos_env', 'gis_handler', 'p_vessel', 'x_disc_pen', 'p_es', 'p_diesel', 'p_pen_power', 'p_feedbarge']

# with wave real data
all_vars_default = all_vars_default  #+ ['p_wave_vec']
# with wave average data
#all_vars_default = all_vars_default  + ['x_wave_ave']

conditions = {'current [m/s]': 'gis/data/Surface Current (m-s).tif',
              'oxygen [mg/l]': 'gis/data/Surface Oxygen (mg-l).tif',
              'salinity [PSU]': 'gis/data/Surface Salinity (PSU).tif',
              'temperature [°C]': 'gis/data/Surface Temperature (C).tif',
              'period [s]': 'gis/data/Wave Energy Period (s).tif',
              'height [m]': 'gis/data/Wave Significant Height (m).tif',
              'bathymetry [m]': 'gis/data/Bathymetry (m).tif',
              #'distance to shore [m]': 'gis/data/Distance to Shore (m).tif',
              'distance to port [m]': 'gis/data/Distance to Port (m).tif'}

# high fishing is above average, very high is more than one standard deviation above average
conflicts = {'very high fishing traffic': 'gis/data/Very High Fishing Vessel Traffic.geojson',
#            'high fishing traffic': 'data/High Fishing Vessel Traffic.geojson',
             'marine protected areas': 'gis/data/Marine Protected Areas.geojson',
             'danger zones': 'gis/data/Danger Zones and Restricted Areas.geojson',
             'submarine': 'gis/data/Submarine Transit Lanes.geojson',
             'torpex': 'gis/data/Cape Cod TORPEX.geojson',
             'block island': 'gis/data/Block Island Renewable Energy Zone.geojson',
             'ma wind': 'gis/data/MA Wind Energy Areas.geojson',
             'wind lease': 'gis/data/Potential Wind Lease Areas.geojson',
             'wind planning': 'gis/data/Wind Planning Areas.geojson',
             'shipping': 'gis/data/Shipping Lanes.geojson'}

waters = "gis/data/Northeast State and Federal Waters.geojson"

handler = GISHandler(conditions, conflicts, waters)

args_default = {}
args_default['wave_data'] = "../Wave Data/32_43.49_-67.88_2009.csv"
args_default['aqua_load'] = "../Aquaculture Load/Load 24 hour.xlsx"
args_default['fixed_num_pen'] = 12
args_default['moo_n_obj'] = 2


def wpaf_single_opt(all_vars_in = None, args_in = None):
    if all_vars_in is not None:
        all_vars = all_vars_in        
    else:
        all_vars = all_vars_default
    
    args = copy.deepcopy(args_default)
    args['handler'] = handler
    if args_in is not None:
        args.update(args_in)

    x_init, p_init, x, soo_res_best, p_best, op_obj, moo_res_best = wpaf_opt(all_vars, args)
    return x_init, p_init, x, soo_res_best, p_best, op_obj


def wpaf_multi_opt(all_vars_in = None, args_in = None):
    if all_vars_in is not None:
        all_vars = all_vars_in        
    else:
        all_vars = all_vars_default
    
    args = copy.deepcopy(args_default)
    args['handler'] = handler
    if args_in is not None:
        args.update(args_in)
    
    args['multi_objective'] = "multi_objective"
    x_init, p_init, x, soo_res_best, p_best, op_obj, moo_res_best = wpaf_opt(all_vars, args)
    return x, moo_res_best, p_best, op_obj


def wpaf_opt(all_vars, args):
    wec_types = ['point absorber'] #,'attenuator','terminator','point_absorber_RM3']

    # design variables and parameters
    if 'p_wave_vec' in all_vars:
        x_name = ['x_wec','x_pen', 'x_es']
        p_name = ['x_type_wec', 'x_env', 'x_wave_env', 'gis_handler', 'pos_env', 'p_pen_power'] #,'x_env'
    else: # es is not required when average wave is given
        x_name = ['x_wec','x_pen']
        p_name = ['x_type_wec', 'x_env', 'x_wave_env', 'x_es', 'gis_handler', 'pos_env', 'p_pen_power'] #,'x_env'

    x = OpData(x_name)
    if 'x0' in args:
        for i in range(len(x.list)):
            x.nom_dict[x.list[i]] = args['x0'][i]

    x_disc_name = ['x_disc_pen']
    x_disc = OpData(x_disc_name)
    if 'fixed_num_pen' in args:
        x_disc.bnds[0] = (args['fixed_num_pen'], args['fixed_num_pen'])   # to accelerate the computation
    
    
    param = OpData(p_name + x_disc_name)

    #Environemental parameter can be set either here or in module.
    #param.nom_dict['temp'] =     16  #'C'
    #param.nom_dict['salinity'] = 33  #'[PSU]'
    #param.nom_dict['U'] =        0.2 #'[m/s]'
    #param.nom_dict['O2_in'] =    8   #'[mg/l]'
    #param.nom_dict['wave_height'] = 1.4  #'[m]'     
    #param.nom_dict['wave_period'] = 8.33 #'[s]'
    #param_val['wave_data'] = "wave_data/Wave_Data2.csv"

    if 'handler' in args:
        param.nom_dict['handler'] = args['handler']
    else:
        print('GIS handler is needed!')
        exit()
    

    if 'wave_data' in args:
        df = pd.read_csv(args['wave_data'])
        param.nom_dict['pos_lat'] = float(df['Latitude'][0])
        param.nom_dict['pos_long'] = float(df['Longitude'][0])

        df = pd.read_csv(args['wave_data'], skiprows=2)
        wave_period = df['Energy Period']
        wave_height = df['Significant Wave Height']
        # noise= np.random.rand(8760)
        # vec = 1.37 * np.ones(8760)
        # vec[1000:1100] = 1.2
        param.nom_dict['wave_height'] = np.array(wave_period.values) #vec + noise/10  #
        param.nom_dict['wave_energy_period'] = np.array(wave_height.values)  #6.44 * np.ones(8760) #

    if 'aqua_load' in args:
        df = pd.read_excel(args['aqua_load'])
        time_index = df['hour']
        time_index = np.array(time_index.values)
        summer_feedbarge_power = df['summer feedbarge kw']
        summer_lighting_power_per_kg = df['summer lighting kw/kg']
        summer_equipment_power_per_kg = df['summer equipment kw/kg']
        winter_feedbarge_power = df['winter feedbarge kw']
        winter_lighting_power_per_kg = df['winter lighting kw/kg']
        winter_equipment_power_per_kg = df['winter equipment kw/kg']

        param.nom_dict['summer_feedbarge_power'] = np.array(summer_feedbarge_power.values)
        param.nom_dict['summer_lighting_power_per_kg'] = np.array(summer_lighting_power_per_kg.values)
        param.nom_dict['summer_equipment_power_per_kg'] = np.array(summer_equipment_power_per_kg.values)
        param.nom_dict['winter_feedbarge_power'] = np.array(winter_feedbarge_power.values)
        param.nom_dict['winter_lighting_power_per_kg'] = np.array(winter_lighting_power_per_kg.values)
        param.nom_dict['winter_equipment_power_per_kg'] = np.array(winter_equipment_power_per_kg.values)

    #optimization
    moo_res = {}
    soo_res = {}
    soo_res_best = {}
    init_flag = 1

    for i in range(len(wec_types)):
        param.nom_dict['wec_type'] = wec_types[i]

        for x_disc_i in range(len(x_disc.list)):
            lower_bnd, upper_bnd = x_disc.bnds[x_disc_i]

            for x_disc_val in range(lower_bnd,upper_bnd+1):
                param.nom_dict[x_disc.list[x_disc_i]] = x_disc_val            

                if 'x_wave_ave' in all_vars:
                    param.nom_dict['es_size'] = 0                    

                if 'multi_objective' not in args:
                    #print('single optimization is running ...')
                    soo_res, op_obj, p = optimization.run_soo_optimization(x.name, x.nom0, param.name, param.nom_dict, all_vars, max_iter)
                else:
                    #print('multi optimization is running ...')
                    moo_res, op_obj, p = optimization.run_moo_optimization(args['moo_n_obj'], x.name, param.name, param.nom_dict, all_vars, max_iter)

                if init_flag:
                    x_init, p_init = x , p
                    soo_res_best = copy.deepcopy(soo_res)
                    moo_res_best = moo_res
                    p_best = p.nom_dict
                elif (soo_res.success) and (soo_res.fun < soo_res_best.fun):
                    soo_res_best = copy.deepcopy(soo_res)
                    p_best = p.nom_dict

                init_flag = 0

    return x_init, p_init, x, soo_res_best, p_best, op_obj, moo_res_best

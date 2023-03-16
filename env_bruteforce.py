from modules import Aqua_Obj
from optimization import OpData
import optimization
from utilities import *
from env_opt import *
import copy
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import time

class GeoData():
    def __init__(self):
        self.points = gpd.GeoDataFrame(columns=['pos_long', 'pos_lat', 'geometry', 'ok-conditions', 'ok-scope', 'ok-conflicts', 'obj_func', 'valid_point'], geometry='geometry')
        
    def store_data(self, valid_point, gis_data, aqua_obj: Aqua_Obj):
        point = Point(float(gis_data["x"]), float(gis_data["y"]))
        data = {'pos_long': float(gis_data["x"]), 'pos_lat': float(gis_data["y"]), 'geometry': point, 
                'ok-conditions': gis_data["ok-conditions"].bool(), 'ok-scope': gis_data["ok-scope"].bool(), 'ok-conflicts': gis_data["ok-conflicts"].bool(),
                'valid_point': gis_data["ok-conditions"].bool() and gis_data["ok-scope"].bool() and gis_data["ok-conflicts"].bool()}
        if valid_point:
            data['obj_func'] = aqua_obj.obj_func
            data['vessel_travel_price [$]'] =  aqua_obj.vessel.price
            data['wec_price [$]'] = aqua_obj.wec.price
            data['fish_yield [kg]'] = aqua_obj.pen.fish_yield
            data['cost_per_yield [$/kg]'] = aqua_obj.cost_per_yield
            data['wec_number [-]'] = aqua_obj.wec.wec_number
            data['aqua_energy [kWh]'] = aqua_obj.pen.power
            data['wec_AEP_per_unit [kWh]'] = aqua_obj.wec.AEP_per_unit
            data['wave_power [kW/m]'] = aqua_obj.wave_in.P_wave
            data['travel_distance [km]'] = aqua_obj.vessel.distance
            data['wec_P_ave [kW]'] = aqua_obj.wec.P_ave
            data['wec_LCOE [$/kWh]'] = aqua_obj.wec.LCOE
            
            const = ['fish_yield_cons', 'env_Umin_cons', 'env_Umax_cons', 'env_tempmin_cons', 'env_tempmax_cons', 'env_salinitymin_cons', 'env_salinitymax_cons', 
                      'env_O2_min_cons', 'env_bathymetry_min_cons', 'env_bathymetry_max_cons']
            for i, const_name in enumerate(const):
                data[const_name]=aqua_obj.ineq_constraint[i]
                data['valid_point'] = data['valid_point'] and (aqua_obj.ineq_constraint[i] >= 0)

        self.points = self.points.append(data, ignore_index=True)
    
    def savefile(self, name):
        self.points.to_file(name + '.geojson', driver='GeoJSON')
        self.points.to_excel(name + '.xlsx') 


def env_bruteforce(all_vars, *args):
    wec_types = ['point absorber'] #,'attenuator','terminator','point_absorber_RM3']

    # design variables
    x_name = ['pos_env']
    x = OpData(x_name)

    # parameters
    param_name = ['x_wec','x_type_wec', 'x_pen', 'gis_handler']
    param = OpData(param_name)

    # WEC and Pen Parameters are defined by optimal results obtained by running "run_sim_wec" 
    param.nom_dict['capture_width']=  np.NaN    #[m]
    param.nom_dict['pen_diameter']=   21.4     #[m]
    param.nom_dict['pen_height']=     19.3     #[m]
    param.nom_dict['stock_density']=  20       #[kg/m^3]

    #print(x.bnds)

    if 'grid_resolution' in args[0]:
        grid_resolution = args[0]['grid_resolution']
    else:
        grid_resolution = 1

    if 'handler' in args[0]:
        handler = args[0]['handler']
        param.nom_dict['handler'] = args[0]['handler']
    else:
        print('GIS handler is needed!')
        exit()

    lat_grid = np.arange(x.bnds[0][0], x.bnds[0][1] + grid_resolution/10, grid_resolution)
    long_grid = np.arange(x.bnds[1][0], x.bnds[1][1] + grid_resolution/10, grid_resolution)

    #brute force
    init_flag_total = 1
    init_flag_without_conflict = 1
    data_history = GeoData()
    best_x_total, aqua_obj_best_total, p_best, best_x_without_conflict, aqua_obj_best_without_conflict= None, None, None, None, None

    for i in range(len(wec_types)):
        param.nom_dict['wec_type'] = wec_types[i]

        # fill default parameters
        p = optimization.argument_fun(x.name, param.name, param.nom_dict, all_vars)

        for pos_lat in (lat_grid):
            for pos_long in (long_grid):
                gis_data = handler.query(pos_long, pos_lat)
                valid_point = gis_data["ok-conditions"].bool() and gis_data["ok-scope"].bool() #and gis_data["ok-conflicts"].bool()
                #print(valid_point, 'lat=', pos_lat, 'long=', pos_long)
                if valid_point:
                    x.nom_dict['pos_lat'] = pos_lat
                    x.nom_dict['pos_long'] = pos_long
                    aqua_obj = Aqua_Obj(x.nom0, x.name, p.nom_dict) 
                    data_history.store_data(valid_point, gis_data, aqua_obj)
                    
                    if sum(ineq_cons < 0 for ineq_cons in aqua_obj.ineq_constraint) == 0:
                        if init_flag_total:
                            best_x_total = copy.deepcopy(x)
                            res_best_total = copy.deepcopy(aqua_obj.obj_func)
                            p_best = copy.copy(p.nom_dict)
                        elif (aqua_obj.obj_func < res_best_total):
                            best_x_total = copy.deepcopy(x)
                            res_best_total = copy.deepcopy(aqua_obj.obj_func)
                            p_best = copy.copy(p.nom_dict)
                        
                        init_flag_total = 0

                        if (gis_data["ok-conflicts"].bool()):
                            if init_flag_without_conflict:
                                best_x_without_conflict = copy.deepcopy(x)
                                res_best_without_conflict = copy.deepcopy(aqua_obj.obj_func)
                            elif (aqua_obj.obj_func < res_best_without_conflict):
                                best_x_without_conflict = copy.deepcopy(x)
                                res_best_without_conflict = copy.deepcopy(aqua_obj.obj_func)
                            init_flag_without_conflict = 0

                    
                else:
                    data_history.store_data(valid_point, gis_data, None)

    timestr = time.strftime("%Y%m%d_%H%M%S")
    filepath = 'results/env_bruteforce_' + timestr + '_' + str(grid_resolution) + '_' + str(aqua_obj.pen.n) + 'cages'
    data_history.savefile(filepath)

    if init_flag_total == 0:
        aqua_obj_best_total = Aqua_Obj(best_x_total.nom0, best_x_total.name, p_best) 
    if init_flag_without_conflict == 0:
        aqua_obj_best_without_conflict = Aqua_Obj(best_x_without_conflict.nom0, best_x_without_conflict.name, p_best) 

    return best_x_total, aqua_obj_best_total, p_best, best_x_without_conflict, aqua_obj_best_without_conflict
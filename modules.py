import numpy as np
import csv
from typing import Tuple
from objects import *
import pandas as pd
from gis.gis_handler import GISHandler

def obj(x_in, x_name, p_in: dict):
    aqua_obj = Aqua_Obj(x_in, x_name, p_in) 
    return aqua_obj.cost_per_yield

        
def ineq_constraint(x_in, x_name, p_in: dict):
    aqua_obj = Aqua_Obj(x_in, x_name, p_in) 
    g = np.array([aqua_obj.P_gen_cons, aqua_obj.fish_yield_cons, 
                  aqua_obj.env_Umin_cons, aqua_obj.env_Umax_cons,
                  aqua_obj.env_tempmin_cons, aqua_obj.env_tempmax_cons, 
                  aqua_obj.env_salinitymin_cons, aqua_obj.env_salinitymax_cons, 
                  aqua_obj.env_O2_min_cons])

    return g


def eq_constraint(x_in, x_name, p_in: dict):
    aqua_obj = Aqua_Obj(x_in, x_name, p_in) 
    h = np.array([0])

    return h
    
class Aqua_Obj(object):
    def __init__(self, x0, x_name, p):
        self.x0 = x0
        self.x_name = x_name
        self.p = p
        
        self.wec, self.pen, self.fish, self.es = input_merge(self.x0, self.x_name, self.p)
                
        self.fish_yield_func()
        self.power()
        self.carrying_capacity = self.pen.carrying_capacity(self.fish)
        self.cost_per_yield = self.price/self.fish_yield 
                
        #self.P_gen_cons = self.wec.annual_energy - self.pen.annual_energy
        self.P_gen_cons = np.min(self.es.P_stored_cum)
        #print('P_gen_cons=', wec.annual_energy, pen.power, P_gen_cons)
        #P_gen_cons = wec.P_gen - pen.power

        # Multiply 0.001 to get a similar order with P_gen_cons
        self.fish_yield_cons = (self.carrying_capacity - self.fish_yield)   
        #fish_yield_cons = (carrying_capacity - fish_yield / pen.n) * 0.001   #for each pen

        self.env_Umin_cons = self.pen.U - self.fish.U_min
        self.env_Umax_cons = self.fish.U_max - self.pen.U
        self.env_tempmin_cons = self.pen.temp - self.fish.temp_min
        self.env_tempmax_cons = self.fish.temp_max - self.pen.temp
        self.env_salinitymin_cons = self.pen.salinity - self.fish.salinity_min
        self.env_salinitymax_cons = self.fish.salinity_max - self.pen.salinity
        self.env_O2_min_cons = self.pen.O2_in - self.fish.O2_min
        #self.es_size_cons = self.es.size - self.es.size_required
    
    def power(self):
        self.es.P_diff = self.wec.P_electrical - self.pen.power_hour
        return
    
    @property
    def price(self):
        price = self.wec.price + self.pen.price + self.fish_feed_price + self.es.price
        return price
    
    #def wave_climate(self):
    #    Hs = self.wec.wave_damping * self.wec.wave_Hs
    #    T = self.wec.wave_T
    #    self.wave_out = Wave(Hs,T)
    #    return
    
    def fish_yield_func(self):
        survival_rate = (1-self.fish.loss_rate) 
        self.fish_yield = self.pen.biomass * survival_rate  # [kg]
        #print("fish_yield", pen.fish_yield)
        return 
    
    @property
    def fish_feed_price(self):    
        return self.pen.biomass * self.fish.FCR * self.fish.feed_unit_cost

    def plot_variable(self):
        self.fish.plot_variable
        return

    def carrying_capacity_print(self):
        return self.pen.TPF_O2, self.carrying_capacity
    
    @property
    def price_breakdown(self):
        return self.wec.price, self.pen.price, self.fish_feed_price
    
    def plot_power(self):
        ax1 = plt.subplot(5,1,1)
        ax1.plot(self.wec.P_mechanical)
        ax1.set(xlabel='time [hour]', ylabel='P_mechanical(kW)');
        ax1.grid()
        #plt.ylim(-10, 10)
        plt.show()

        ax1 = plt.subplot(5,1,2)
        ax1.plot(self.wec.P_electrical)
        ax1.set(xlabel='time [hour]', ylabel='P_electrical(kW)');
        ax1.grid()
        #plt.ylim(0, 20)
        plt.show()

        ax1 = plt.subplot(5,1,3)
        ax1.plot(self.pen.power_hour)
        ax1.set(xlabel='time [hour]', ylabel='pen power(kW)');
        ax1.grid()
        #plt.ylim(-10, 10)
        plt.show()

        ax1 = plt.subplot(5,1,4)
        ax1.plot(self.es.P_diff)
        ax1.set(xlabel='time [hour]', ylabel='P_diff(kW)');
        ax1.grid()
        #plt.ylim(-20, 20)
        plt.show()

        '''
        ax2 = plt.subplot(5,1,4)
        ax2.plot(self.es.P_stored_hour, label='kWh')
        ax2.set(xlabel='time [hour]', ylabel='P_stored_hour(kWh)');
        ax2.legend()
        plt.show()
        '''

        ax3 = plt.subplot(5,1,5)
        ax3.plot(self.es.P_stored_cum)
        ax3.set(xlabel='time [hour]', ylabel='P_stored_cum (kW)');
        ax3.grid()
        plt.show()


def input_merge(x_in, x_name, p):
    # merge input dicts
    
    x_list = variable_lookup(x_name)
    
    if 'x0_scale' in p:
        scale = p['x0_scale']
    else:
        scale = np.ones(len(x_list))
    
    x = {}
    for i in range(len(x_list)):
        x[x_list[i]] = x_in[i]*scale[i]
 
    if 'pos_env' in x_name:
        gis_data = import_gis_data(x_in[0], x_in[1])
        p['U'] = float(gis_data["current"])
        p['O2_in'] = float(gis_data["oxygen"])
        p['salinity'] = float(gis_data["salinity"])
        p['temp'] = float(gis_data["temperature"])

    if ('wave_data' in p) and (p['wave_data']!=""):
        wave_file_name = p['wave_data']
        wave_period_i, wave_height_i = load_wave_data(wave_file_name)
        p['wave_height'] = wave_height_i
        p['wave_period'] = wave_period_i
    elif 'wave_height' in p:
        p['wave_height'] = p['wave_height'] * np.ones(8760)
        p['wave_period'] = p['wave_period'] * np.ones(8760)
    elif 'wave_height' in x:
        x['wave_height'] = x['wave_height'] * np.ones(8760)
        x['wave_period'] = x['wave_period'] * np.ones(8760)

    ins = {**x, **p}

    #print(x)
    #print(p)
    # create objects    
    wec = WEC(ins['capture_width'], ins['capture_width_ratio_dict'],
            ins['wave_damping_dict'], ins['wec_type'], ins['wec_unit_cost'],
            ins['wave_height'], ins['wave_period'],
            ins['eta'], ins['capacity_factor'])

    pen = Pen(ins['pen_diameter'], ins['pen_height'], ins['pen_depth'], ins['stock_density'], 
              ins['num_pens'], ins['spacing'], ins['pen_unit_cost'], ins['temp'], 
              ins['O2_in'], ins['U'], ins['salinity'], ins['permeability'])
    
    fish = Fish(ins['F_f'], ins['F_p'], ins['F_c'], ins['A_f'], ins['A_p'], ins['A_c'],
                ins['O_f'], ins['O_p'], ins['O_c'], ins['C_f'], ins['C_p'], ins['C_c'],
                ins['P_f'], ins['P_p'], ins['tau'], ins['loss_rate'], ins['harvest_weight'], 
                ins['O2_min'], ins['U_min'], ins['U_max'], ins['temp_min'], ins['temp_max'], 
                ins['salinity_min'], ins['salinity_max'], ins['FCR'], ins['feed_unit_cost'])
    
    es = ES(ins['es_eta'], ins['es_dod'], ins['es_unit_cost'], ins['es_size'])

    #print(ins['es_size'])
    return wec, pen, fish, es

def variable_lookup(var_category_names):
    # input is a cell array containing some subset of the following categories:
    # x_wec, x_wec_type, x_pen, p_pen, p_env, p_wec, p_fish
    # output is a cell array with variable names matching the input categoeies.

    var_list = []
    if any('x_wec' in i for i in var_category_names):
        var_list.append('capture_width')

    if any('x_type_wec' in i for i in var_category_names):
        var_list.append('wec_type')
        
    if any('x_pen' in i for i in var_category_names):
        var_list.append('pen_diameter')
        var_list.append('pen_height')
        var_list.append('spacing')
        var_list.append('stock_density')
        var_list.append('pen_depth')
    
    if any('x_disc_pen' in i for i in var_category_names):
        var_list.append('num_pens')

    if any('p_pen' in i for i in var_category_names):
        var_list.append('pen_unit_cost')
        var_list.append('permeability')
        
    if any('pos_env' in i for i in var_category_names):
        var_list.append('pos_x')
        var_list.append('pos_y')

    if any('x_env' in i for i in var_category_names):
        var_list.append('temp')
        var_list.append('O2_in')
        var_list.append('salinity')
        var_list.append('U')
        var_list.append('wave_height')
        var_list.append('wave_period')
    
    if any('wave_data' in i for i in var_category_names):
        var_list.append('wave_data')

    if any('p_wec' in i for i in var_category_names):
        var_list.append('wec_unit_cost')
        var_list.append('capture_width_ratio_dict')
        var_list.append('wave_damping_dict')
        var_list.append('eta')
        var_list.append('capacity_factor')
    
    if any('p_fish' in i for i in var_category_names):
        var_list.append('F_f')             #Fraction of Fat in the Feed Mix [-]
        var_list.append('F_p')             #Fraction of Protein in the Feed Mix [-]
        var_list.append('F_c')             #Fraction of Carbohydrates in the Feed Mix [-]
        var_list.append('A_f')             #Assimilated Fraction of Fat [-]
        var_list.append('A_p')             #Assimilated Fraction of Protein [-]
        var_list.append('A_c')             #Assimilated Fraction of Carbohydrate [-]
        var_list.append('O_f')             #Oxygen Demand to Break Down Fat [gO2/gFat]
        var_list.append('O_p')             #Oxygen Demand to Break Down Protein [gO2/gProtein]
        var_list.append('O_c')             #Oxygen Demand to Break Down Carbohydrate [gO2/gCarbohydrate]
        var_list.append('C_f')             #Fat's Specific Energy Content [cal/g]
        var_list.append('C_p')             #Protein's Specific Energy Content [cal/g]
        var_list.append('C_c')             #Carbohydrate's Specific Energy Content [cal/g]
        var_list.append('P_f')             #Fat Content of Fish [-]
        var_list.append('P_p')             #Protein Content of Fish [-]
        var_list.append('tau')             #Inverse Temperature Scale [1/C]
        var_list.append('loss_rate')       #Fish Loss Rate [-]
        var_list.append('harvest_weight')  #Fish Harvest Size [kg/fish]
        var_list.append('O2_min')          #Dissolved Oxygen Threshold [%]
        var_list.append('U_min')
        var_list.append('U_max')
        var_list.append('temp_min')
        var_list.append('temp_max')
        var_list.append('salinity_min')
        var_list.append('salinity_max')
        var_list.append('FCR')             #Feed Conversion Ratio [kgFeed/kgFish]
        var_list.append('feed_unit_cost')  #[$/kgFeed]
    
    if any('x_es' in i for i in var_category_names):
        var_list.append('es_size')

    if any('p_es' in i for i in var_category_names):
        var_list.append('es_eta')
        var_list.append('es_dod')
        var_list.append('es_unit_cost')

    if len(var_list)==0:
        print('Your input did not match any of the category names.', var_category_names)
    

    return var_list


def default_values(var_category_names):
    vals = {}
    wec_types = (['attenuator','terminator','point_absorber_RM3'], '[-]')
    capture_width_ratios = ([0.16, 0.34, 0.16], '[-]') 
    wave_dampings = ([0, 0.13, 0.17], '[-]')      

    if any('x_wec' in i for i in var_category_names):
        vals['capture_width'] = (4, '[m]')     #12

    if any('x_type_wec' in i for i in var_category_names):
        vals['wec_type'] = ('point_absorber_RM3', '[-]')
        
    if any('x_pen' in i for i in var_category_names):
        vals['pen_diameter'] = (15, '[m]')    #20 
        vals['pen_height'] = (4, '[m]')      #6  
        vals['spacing'] = (150, '[m]')          
        vals['stock_density'] = (20 , '[kg/m^3]') #20
        vals['pen_depth'] = (10, '[m]')         
   
    if any('x_disc_pen' in i for i in var_category_names):
        vals['num_pens'] = (18, '[-]')          
        
    if any('p_pen' in i for i in var_category_names):
        vals['pen_unit_cost'] = (100, '[$/m^3]')    # 80 $/m^3 for net pen + 20 $/m^3 for mooring
        vals['permeability'] = (0.8, '[-]')
        
    if any('pos_env' in i for i in var_category_names):
        vals['pos_x'] = (-70, '[m]')
        vals['pos_y'] = (41, '[m]')
    
    if any('x_env' in i for i in var_category_names):
        vals['temp'] = (16, 'C')               
        vals['salinity'] = (33, '[PSU]')
        vals['U'] = (.2, '[m/s]')
        vals['O2_in'] = (8,'[mg/l]')
        vals['wave_height'] = (1.4, '[m]')
        vals['wave_period'] = (8.33, '[s]')
    
    if any('wave_data' in i for i in var_category_names):
        vals['wave_data'] = ("",'[-]')
    
    if any('p_wec' in i for i in var_category_names):
        vals['wec_unit_cost'] = (0.45 * 1.19, '[$/kWh]')   # 'point_absorber_RM3' * inflation rate from 2014 to 2022
        vals['capture_width_ratio_dict'] = (dict(zip(wec_types[0], capture_width_ratios[0])), '[-]')
        vals['wave_damping_dict'] = (dict(zip(wec_types[0], wave_dampings[0])), '[-]')
        vals['eta'] = (0.8,'[-]') 
        vals['capacity_factor'] = (0.3,'[-]') 
        
    if any('p_fish_salmon' in i for i in var_category_names):
        vals['F_p'] = (0.45, '[-]')
        vals['F_c'] = (0.07, '[-]')
        vals['F_f'] = (0.15, '[-]')
        vals['A_p'] = (0.97, '[-]') #0.89
        vals['A_c'] = (0.6, '[-]') #0.5
        vals['A_f'] = (0.9, '[-]')
        vals['O_p'] = (1.89, '[gO2/gProtein]')
        vals['O_c'] = (1.07, '[gO2/gCarbohydrate]')
        vals['O_f'] = (2.91, '[gO2/gFat]')
        vals['C_p'] = (5650, '[cal/g]')
        vals['C_c'] = (4100, '[cal/g]')
        vals['C_f'] = (9450, '[cal/g]')
        vals['P_f'] = (0.18, '[-]')
        vals['P_p'] = (0.18, '[-]')
        vals['tau'] = (0.08, '[1/C]')
        vals['loss_rate'] = (0.15, '[-]')
        vals['harvest_weight'] = (4, '[kg/fish]') # can be between 4 [kg] to 6 [kg] 
        vals['O2_min'] = (4.41, '[mg/l]')
        vals['U_min'] = (0.01,'[m/s]') # From Environment
        vals['U_max'] = (2,'[m/s]')
        vals['temp_min'] = (2,'[C]')
        vals['temp_max'] = (28,'[C]')
        vals['salinity_min'] = (30,'[PSU]')
        vals['salinity_max'] = (35,'[PSU]')
        vals['FCR'] = (1.15,'[kgFeed/kgFish]')
        vals['feed_unit_cost'] = (1.48,'[$/kgFeed]')
        
    if any('x_es' in i for i in var_category_names):
        vals['es_size'] = (50,'[kWh]')
    
    if any('p_es' in i for i in var_category_names):
        vals['es_eta'] = (0.92,'[-]')
        vals['es_dod'] = (0.7,'[-]')
        vals['es_unit_cost'] = (271,'[$/kWh]')

    
    '''
    if any('p_fish_black_sea_bass' in i for i in var_category_names):
        vals['F_p'] = 2
        vals['F_c'] = 2
        vals['F_f'] = 2
        vals['A_p'] = 2
        vals['A_c'] = 2
        vals['A_f'] = 2
        vals['O_p'] = 2
        vals['O_c'] = 2
        vals['O_f'] = 2
        vals['C_p'] = 2000
        vals['C_c'] = 2000
        vals['C_f'] = 20000
        vals['P_f'] = 2
        vals['P_p'] = 2
        vals['O2_min'] = 2
        vals['tau'] = 2
        vals['loss_rate'] = 2
        vals['harvest_weight'] = 2
    '''
    
    #assert(fieldnames(vals) == variable_lookup(var_category_names));
    return vals

def bnds_values(var_category_names):
    bnds = {}

    if any('x_wec' in i for i in var_category_names):
        bnds['capture_width'] = (1, 40)     #[m] 40
    
    if any('x_pen' in i for i in var_category_names):
        bnds['pen_diameter'] = (10, 40)      #[m]
        bnds['pen_height'] = (3, 30)        #[m] (5,30)
        bnds['spacing'] = (100, 300)         #[m]
        bnds['stock_density'] = (15, 30)     #[kg/m^3]
        bnds['pen_depth'] = (1, 30)         #[m] 

    if any('x_disc_pen' in i for i in var_category_names):
        bnds['num_pens'] = (1, 20)         #[-] 
        
    if any('pos_env' in i for i in var_category_names):
        bnds['pos_x'] = (-100, 100)         #[m]
        bnds['pos_y'] = (39, 100)         #[m]
        
    if any('x_env' in i for i in var_category_names):
        bnds['temp'] = (7, 22)              #[C]
        bnds['salinity'] = (26, 36)         #[PSU]
        bnds['O2_in'] = (7, 10)             #[mg/l]
        bnds['U'] = (0.01, 0.55)            #[m/s]
        bnds['wave_height'] = (0.2, 3)      #[m]
        bnds['wave_period'] = (1, 12)       #[s]
    
    if any('x_es' in i for i in var_category_names):
        bnds['es_size'] = (5, 500)          #[kWh]
    
    return bnds

def import_gis_data(pos_x, pos_y):
    # Import raster file
    raster_files = {'current': 'gis/data/surface-current-ms.tif',
                'oxygen': 'gis/data/surface-oxygen-mgpl.tif',
                'salinity': 'gis/data/surface-salinity-psu.tif',
                'temperature': 'gis/data/surface-temperature-c.tif'}
    handler = GISHandler(raster_files)
    return handler.query(pos_x, pos_y)

def load_wave_data(file_name):
    df = pd.read_csv(file_name)
    wave_period = df['Peak Period']
    wave_height = df['Significant Wave Height']
    return np.array(wave_period.values), np.array(wave_height.values)
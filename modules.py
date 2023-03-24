import numpy as np
from objects import * 
import math 

def obj(x_in, x_name, p_in: dict):
    aqua_obj = Aqua_Obj(x_in, x_name, p_in) 
    if aqua_obj.valid_point:
        return aqua_obj.obj_func #aqua_obj.cost_per_yield
    else:
        return aqua_obj.obj_func + 30000

def ineq_constraint(x_in, x_name, p_in: dict):
    aqua_obj = Aqua_Obj(x_in, x_name, p_in) 
    if aqua_obj.valid_point: 
        g = np.array(aqua_obj.ineq_constraint)
    else: 
        g = np.array([-1])
    return g


def eq_constraint(x_in, x_name, p_in: dict):
    aqua_obj = Aqua_Obj(x_in, x_name, p_in) 
    if aqua_obj.valid_point: 
        h = np.array([0])
    else:
        h = np.array([-1])
    return h

def obj_terms(x_in, x_name, p_in: dict):
    aqua_obj = Aqua_Obj(x_in, x_name, p_in) 
    if aqua_obj.valid_point: 
        return np.array([aqua_obj.cost_NPV, aqua_obj.pen.fish_yield, aqua_obj.pen.price, aqua_obj.wec.price])
    else:
        return np.array([-1, -1, -1, -1])

class Aqua_Obj(object):
    def __init__(self, x0, x_name, p):
        self.x0 = x0
        self.x_name = x_name
        self.p = p
        
        self.wec, self.wave_in, self.pen, self.fish, self.vessel, self.valid_point, self.gis_data = input_merge(self.x0, self.x_name, self.p)
        
        if self.valid_point:
            self.fish_yield_func()

            if math.isnan(self.wec.capture_width):
                self.wec.capture_width = self.wec.set_capture_width(self.pen, self.wave_in)

            self.wec.P_gen = self.power()

            #print('    ', self.vessel.price, self.wec.price, self.wec.LCOE, self.wec.LCOE_base_RM3, self.wec.AEP_per_unit, self.wec.AEP, self.wec.capture_width, self.wec.wec_number,self.pen.power, self.wave_in.P_wave)

            self.carrying_capacity = self.pen.carrying_capacity(self.fish)
            self.cost_per_yield = self.cost_NPV/self.pen.fish_yield 

            # power supply constraint to ensure supply power demand of net pen
            #self.P_gen_cons = (self.wec.AEP - self.pen.power) / self.wec.AEP

            # fish yield constraint to ensure a healthy offshore environment
            self.fish_yield_cons = self.carrying_capacity - self.pen.fish_yield
        
            # net pen geometry constraint to present a practical design and ratio between diameter and height
            #self.pen_ratio_low_cons = (self.pen.D - self.pen.H) / self.pen.D
            #self.pen_ratio_up_cons = (3*self.pen.H - self.pen.D) / (3*self.pen.H)

            self.env_Umin_cons = self.pen.U - self.fish.U_min
            self.env_Umax_cons = self.fish.U_max - self.pen.U
            self.env_tempmin_cons = self.pen.temp - self.fish.temp_min
            self.env_tempmax_cons = self.fish.temp_max - self.pen.temp
            self.env_salinitymin_cons = self.pen.salinity - self.fish.salinity_min
            self.env_salinitymax_cons = self.fish.salinity_max - self.pen.salinity
            self.env_O2_min_cons = self.pen.O2_in - self.fish.O2_min
            self.env_bathymetry_min_cons = self.pen.bathymetry - self.pen.H - self.pen.waterdepth_underpen_min
            self.env_bathymetry_max_cons = self.pen.waterdepth_underpen_max + self.pen.H - self.pen.bathymetry
        
    @property
    def obj_func(self):
        return self.cost_NPV
    
    @property
    def ineq_constraint(self):
        return [self.fish_yield_cons , self.env_Umin_cons, self.env_Umax_cons,
                self.env_tempmin_cons, self.env_tempmax_cons, 
                self.env_salinitymin_cons, self.env_salinitymax_cons, 
                self.env_O2_min_cons, self.env_bathymetry_min_cons, self.env_bathymetry_max_cons]

    @property
    def cost_NPV(self): #net present value
        cost_NPV = self.wec.cost_NPV + self.vessel.cost_NPV
        return cost_NPV
    
   # @property
   # def price(self):
   #     price = self.wec.price + self.pen.price + self.vessel.price + self.fish_feed_price
   #     return price
    
    def fish_yield_func(self):
        survival_rate = (1-self.fish.loss_rate) 
        self.pen.fish_yield = self.pen.biomass * survival_rate 
        return 
    
    def power(self):
        return self.wec.P_electrical(self.wave_in)
    
    @property
    def fish_feed_price(self):    
        return self.pen.biomass * self.fish.FCR * self.fish.feed_unit_cost
    
    def plot_variable(self):
        self.fish.plot_variable
        return
    
    def carrying_capacity_print(self):
        return self.pen.TPF_O2, self.carrying_capacity

def input_merge(x_in, x_name, p):
    # merge input dicts
    
    x_list = variable_lookup(x_name)
    x = {}
    for i in range(len(x_list)):
        x[x_list[i]] = x_in[i]
 
    gis_data = []
    if 'pos_env' in x_name:
        if 'handler' in p:
            #print('lat=', x_in[0], 'long=', x_in[1], end='              ')
            gis_data =  p['handler'].query(x_in[1], x_in[0]) #import_gis_data(x_in[0], x_in[1])
        else:
            print('GIS handler is needed!')
            exit()
        p['U'] = float(gis_data["current [m/s]"])
        p['O2_in'] = float(gis_data["oxygen [mg/l]"])
        p['salinity'] = float(gis_data["salinity [PSU]"])
        p['temp'] = float(gis_data["temperature [°C]"])
        p['bathymetry'] = -float(gis_data["bathymetry [m]"])
        p['wave_energy_period'] = float(gis_data["period [s]"])
        p['wave_height'] = float(gis_data["height [m]"])
        p['distance'] = float(gis_data["distance to port [m]"]) / 1000
        #p[] = float(gis_data["distance to shore [m]"])
        valid_point = gis_data["ok-conditions"].bool() and gis_data["ok-scope"].bool() #and gis_data["ok-conflicts"].bool()
        #print(valid_point)
    else:
        valid_point = True

    #print(p['wave_height'], p['wave_energy_period'])
    ins = {**x, **p}

    # create objects
    wave_in = Wave(ins['wave_height'], ins['wave_energy_period'])
    
    wec = WEC(ins['capture_width'], ins['capture_width_ratio_dict'],
            ins['wave_damping_dict'], ins['wec_type'],
            ins['capacity_factor'], ins['eta'], ins['float_diameter'],
            ins['CapEx_ref'], ins['OpEx_ref'], ins['lifetime'], ins['discount_rate'])

    pen = Pen(ins['pen_diameter'], ins['pen_height'], ins['pen_depth'], ins['stock_density'], 
              ins['num_pens'], ins['spacing'], ins['pen_unit_cost'], ins['temp'], 
              ins['O2_in'], ins['U'], ins['salinity'], ins['permeability'], ins['bathymetry'],
              ins['pos_lat'], ins['pos_long'])
    
    fish = Fish(ins['F_f'], ins['F_p'], ins['F_c'], ins['A_f'], ins['A_p'], ins['A_c'],
                ins['O_f'], ins['O_p'], ins['O_c'], ins['C_f'], ins['C_p'], ins['C_c'],
                ins['P_f'], ins['P_p'], ins['tau'], ins['loss_rate'], ins['harvest_weight'], 
                ins['O2_min'], ins['U_min'], ins['U_max'], ins['temp_min'], ins['temp_max'], 
                ins['salinity_min'], ins['salinity_max'], ins['FCR'], ins['feed_unit_cost'])
    
    vessel = Vessel(ins['vessel_fuel_consump_rate'], ins['vessel_fuel_cost'], 
                ins['captain_salary'], ins['crew_salary'], ins['crew_num'], 
                ins['time_feed'], ins['vessel_velocity'], ins['travel_number'], ins['distance'],
                ins['lifetime'], ins['discount_rate'])
                
    return wec, wave_in, pen, fish, vessel, valid_point, gis_data


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
        var_list.append('stock_density')

    if any('p_pen' in i for i in var_category_names):
        var_list.append('num_pens')
        var_list.append('spacing')
        var_list.append('pen_depth')
        var_list.append('pen_unit_cost')
        var_list.append('permeability')
        
    if any('pos_env' in i for i in var_category_names):
        var_list.append('pos_lat')
        var_list.append('pos_long')
    
    if any('gis_handler' in i for i in var_category_names):
        var_list.append('handler')
        
    if any('x_env' in i for i in var_category_names):
        var_list.append('temp')
        var_list.append('O2_in')
        var_list.append('wave_height')
        var_list.append('wave_energy_period')
        
    if any('p_env' in i for i in var_category_names):
        var_list.append('salinity')
        var_list.append('U')
        var_list.append('bathymetry')  # water_depth
        
    if any('p_wec' in i for i in var_category_names):
        var_list.append('capture_width_ratio_dict')
        var_list.append('wave_damping_dict')
        var_list.append('eta')
        var_list.append('capacity_factor')
        var_list.append('float_diameter')
        var_list.append('CapEx_ref')
        var_list.append('OpEx_ref')
        var_list.append('lifetime')
        var_list.append('discount_rate')
    
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
    
    if any('p_vessel' in i for i in var_category_names):
        var_list.append('vessel_velocity')
        var_list.append('vessel_fuel_consump_rate')
        var_list.append('vessel_fuel_cost')
        var_list.append('captain_salary')
        var_list.append('crew_salary')
        var_list.append('crew_num')
        var_list.append('time_feed')
        var_list.append('travel_number')
        var_list.append('distance')

    if len(var_list)==0:
        print('Your input did not match any of the category names.', var_category_names)
    

    return var_list


def default_values(var_category_names):
    vals = {}
    wec_types = (['attenuator','terminator','point absorber'], '[-]')
    capture_width_ratios = ([0.16, 0.34, 0.16], '[-]') 
    wave_dampings = ([0, 0.13, 0.17], '[-]')            

    if any('x_wec' in i for i in var_category_names):
        vals['capture_width'] = (20, '[m]')

    if any('x_type_wec' in i for i in var_category_names):
        vals['wec_type'] = ('point absorber', '[-]')
        
    if any('x_pen' in i for i in var_category_names):
        vals['pen_diameter'] = (21.4, '[m]') 
        vals['pen_height'] = (19.3, '[m]')  
        vals['stock_density'] = (20 , '[kg/m^3]')
   
    if any('p_pen' in i for i in var_category_names):
        vals['num_pens'] = (12, '[-]')  #{5, 12, 40}
        vals['spacing'] = (150, '[m]')
        vals['pen_depth'] = (10, '[m]')   
        vals['pen_unit_cost'] = (100, '[$/m^3]')    # 80 $/m^3 for net pen + 20 $/m^3 for mooring
        vals['permeability'] = (0.8, '[-]')      
        
    if any('pos_env' in i for i in var_category_names):
        vals['pos_lat'] = (42.0, 'm')  #42.203
        vals['pos_long'] = (-70.0, 'm') #70.154
    
    if any('gis_handler' in i for i in var_category_names):
        vals['handler'] = ({}, '[-]') 

    if any('x_env' in i for i in var_category_names):
        vals['temp'] = (10.29, 'C') 
        vals['O2_in'] = (9.5,'[mg/l]')
        vals['wave_height'] = (1.37, '[m]')  
        vals['wave_energy_period'] = (6.44, '[s]')
    
    if any('p_env' in i for i in var_category_names):
        vals['salinity'] = (31.6, '[PSU]')
        vals['U'] = (.1, '[m/s]')
        vals['bathymetry'] = [-16.2, '[m]']
    
    if any('p_wec' in i for i in var_category_names):
        vals['capture_width_ratio_dict'] = (dict(zip(wec_types[0], capture_width_ratios[0])), '[-]')
        vals['wave_damping_dict'] = (dict(zip(wec_types[0], wave_dampings[0])), '[-]')
        vals['eta'] = (0.8,'[-]') 
        vals['capacity_factor'] = (0.3,'[-]') 
        vals['float_diameter'] = (20, '[m]')
        vals['CapEx_ref'] = (4833448, '[$]')
        vals['OpEx_ref'] = (116050, '[$]')
        vals['lifetime'] = (20, '[year]')
        vals['discount_rate'] = (0.07, '[-]')
        
    if any('p_fish_salmon' in i for i in var_category_names):
        vals['F_f'] = (0.15, '[-]')                     #Fraction of Fat in the Feed Mix [-]
        vals['F_p'] = (0.45, '[-]')                     #Fraction of Protein in the Feed Mix [-]
        vals['F_c'] = (0.07, '[-]')                     #Fraction of Carbohydrates in the Feed Mix [-]
        vals['A_f'] = (0.9, '[-]')                      #Assimilated Fraction of Fat [-]
        vals['A_p'] = (0.97, '[-]')                     #Assimilated Fraction of Protein [-]
        vals['A_c'] = (0.6, '[-]')                      #Assimilated Fraction of Carbohydrate [-]
        vals['O_f'] = (2.91, '[gO2/gFat]')              #Oxygen Demand to Break Down Fat [gO2/gFat]
        vals['O_p'] = (1.89, '[gO2/gProtein]')          #Oxygen Demand to Break Down Protein [gO2/gProtein]
        vals['O_c'] = (1.07, '[gO2/gCarbohydrate]')     #Oxygen Demand to Break Down Carbohydrate [gO2/gCarbohydrate]
        vals['C_f'] = (9450, '[cal/g]')                 #Fat's Specific Energy Content [cal/g]
        vals['C_p'] = (5650, '[cal/g]')                 #Protein's Specific Energy Content [cal/g]
        vals['C_c'] = (4100, '[cal/g]')                 #Carbohydrate's Specific Energy Content [cal/g]
        vals['P_f'] = (0.18, '[-]')                     #Fat Content of Fish [-]
        vals['P_p'] = (0.18, '[-]')                     #Protein Content of Fish [-]
        vals['tau'] = (0.08, '[1/C]')                   #Inverse Temperature Scale [1/C]
        vals['loss_rate'] = (0.15, '[-]')               #Fish Loss Rate [-]
        vals['harvest_weight'] = (4, '[kg/fish]')       #Fish Harvest Size [kg/fish]  # can be between 4 [kg] to 6 [kg] 
        vals['O2_min'] = (4.41, '[mg/l]')               #Dissolved Oxygen Threshold [%]
        vals['U_min'] = (0.01,'[m/s]')                  # From NorthEast U.S. Environment
        vals['U_max'] = (2,'[m/s]')
        vals['temp_min'] = (2,'[C]')
        vals['temp_max'] = (20,'[C]')
        vals['salinity_min'] = (30,'[PSU]')
        vals['salinity_max'] = (35,'[PSU]')
        vals['FCR'] = (1.5,'[kgFeed/kgFish]')          #Feed Conversion Ratio [kgFeed/kgFish]
        vals['feed_unit_cost'] = (1.48,'[$/kgFeed]')   #[$/kgFeed]
    
    if any('p_vessel' in i for i in var_category_names):
        vals['vessel_velocity'] = (10 * 1.852,'[km/h]')  # 10 knots
        vals['vessel_fuel_consump_rate'] = (25,'[gal/h]')
        vals['vessel_fuel_cost'] = (4.3,'[$/gal]')
        vals['captain_salary'] = (57,'[$/h]')
        vals['crew_salary'] = (32,'[$/h]')
        vals['crew_num'] = (2,'[-]')
        vals['time_feed'] = (2,'[h]')
        vals['travel_number'] = (52,'[-]') # once per week for a year
        vals['distance'] = (76157.726562/1000, '[km]')
    
    #assert(fieldnames(vals) == variable_lookup(var_category_names));
    return vals

def bnds_values(var_category_names):
    bnds = {}

    if any('x_wec' in i for i in var_category_names):
        bnds['capture_width'] = (1, 23)      #[m]  
    
    if any('x_pen' in i for i in var_category_names):
        bnds['pen_diameter'] = (10, 45)      #[m]  
        bnds['pen_height'] = (10, 30)        #[m]  
        bnds['stock_density'] = (10, 20)     #[kg/m^3]
    
    if any('pos_env' in i for i in var_category_names):
        bnds['pos_lat'] = (38.4, 45.2)        #[m]
        bnds['pos_long'] = (-75.8, -65.7)         #[m]
    
    if any('x_env' in i for i in var_category_names):
        bnds['temp'] = (1, 50)              #[C] 
        bnds['O2_in'] = (1, 50)             #[mg/l]
        bnds['wave_height'] = (0.2, 3)      #[m]
        bnds['wave_period'] = (1, 12)       #[s]
    
    return bnds
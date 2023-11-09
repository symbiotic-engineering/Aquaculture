import numpy as np
from objects import * 
import math 

def obj(x_in, x_name, p_in):
    wpaf = WPAF(x_in, x_name, p_in) 
    if wpaf.valid_point:
        return wpaf.obj_func 
    else:
        return wpaf.obj_func + 30000

# multi-objective function
def multi_obj(x_in, x_name, p_in: dict):
    wpaf = WPAF(x_in, x_name, p_in) 
    if wpaf.valid_point:
        return np.array(wpaf.multi_obj_func)
    else:
        return np.array([-1, -1])
    
def ineq_constraint(x_in, x_name, p_in):
    wpaf = WPAF(x_in, x_name, p_in) 
    if wpaf.valid_point: 
        g = np.array(wpaf.ineq_constraint)
    else: 
        g = np.array([-1])
    return g


def eq_constraint(x_in, x_name, p_in):
    wpaf = WPAF(x_in, x_name, p_in) 
    if wpaf.valid_point: 
        h = np.array([0])
    else:
        h = np.array([-1])
    return h

def obj_terms(x_in, x_name, p_in):
    wpaf = WPAF(x_in, x_name, p_in) 
    if wpaf.valid_point: 
        return np.array([wpaf.obj_func, wpaf.cost_per_yield, wpaf.cost_NPV, wpaf.aqua.fish_yield, wpaf.es.total_size]) #, wpaf.aqua.cost_NPV, wpaf.wec.price
    else:
        return np.array([-1, -1, -1, -1, -1, -1])

# ============================================================================ #
#                       WPAF Problem Definition (OF, Constraint,...)           #
# ============================================================================ #

class WPAF(object):
    def __init__(self, x0, x_name, p):
        self.x0 = x0
        self.x_name = x_name
        self.p = p
        
        self.wec, self.aqua, self.es, self.dieselgen, self.valid_point, self.gis_data = input_merge(self.x0, self.x_name, self.p)
        
        if self.valid_point:
            #self.fish_yield_func()
            #print(self.aqua.power)

            # if math.isnan(self.wec.capture_width):
            #     self.wec.capture_width = self.wec.set_capture_width(self.pen, self.wave_in)

            self.wec.P_gen = self.power()

            self.es.sizing_func(self.wec.P_gen - self.aqua.power)

            #self.carrying_capacity = self.aqua.carrying_capacity(self.fish)

            # power supply constraint to ensure supply power demand of net pen
            self.P_gen_cons = 1# min((self.wec.P_gen + self.es.power - self.aqua.power) / self.wec.P_gen)

            # fish yield constraint to ensure a healthy offshore environment
            self.fish_yield_cons = (self.aqua.carrying_capacity - self.aqua.fish_yield) / self.aqua.carrying_capacity

            self.env_Umin_cons = 1 #self.aqua.U - self.fish.U_min
            self.env_Umax_cons = 1 #self.fish.U_max - self.aqua.U
            self.env_tempmin_cons = 1 #self.aqua.temp - self.fish.temp_min
            self.env_tempmax_cons = 1 #self.fish.temp_max - self.aqua.temp
            self.env_salinitymin_cons = 1 #self.aqua.salinity - self.fish.salinity_min
            self.env_salinitymax_cons = 1 #self.fish.salinity_max - self.aqua.salinity
            self.env_O2_min_cons = 1 #self.aqua.O2_in - self.fish.O2_min
            self.env_bathymetry_min_cons = 1 #self.aqua.bathymetry - self.aqua.H - self.aqua.waterdepth_underpen_min
            self.env_bathymetry_max_cons = 1 #self.aqua.waterdepth_underpen_max + self.aqua.H - self.aqua.bathymetry

            # self.pen_ratio_low_cons = (self.aqua.netpen.D - 1.5*self.aqua.netpen.H) / self.aqua.netpen.D
            # self.pen_ratio_up_cons = (2*self.aqua.netpen.H - self.aqua.netpen.D) / (2*self.aqua.netpen.H)
            # self.pen_ratio_cons = 1 #(self.aqua.netpen_geometry(self.aqua.netpen.D) - self.aqua.netpen.H) / self.aqua.netpen.D

            self.sustainable_power_operation_cons = (np.mean(self.wec.P_gen) - np.mean(self.aqua.power)) / np.mean(self.wec.P_gen)
        
    @property
    def obj_func(self):
        return self.cost_per_yield # -self.aqua.fish_yield / 1000000 #self.cost_NPV/100000000 #cost_per_yield
    
    @property
    def multi_obj_func(self):
        #return self.cost_per_yield, -self.aqua.fish_yield / 1000000
        return self.cost_NPV / 100000000, -self.aqua.fish_yield / 1000000
    
    @property
    def ineq_constraint(self):
        return [#self.P_gen_cons, 
                self.fish_yield_cons, 
                # self.env_Umin_cons, self.env_Umax_cons,
                # self.env_tempmin_cons, self.env_tempmax_cons, 
                # self.env_salinitymin_cons, self.env_salinitymax_cons, 
                # self.env_O2_min_cons, self.env_bathymetry_min_cons, self.env_bathymetry_max_cons,
                # self.pen_ratio_low_cons, self.pen_ratio_up_cons,
                self.sustainable_power_operation_cons]
    
    # @property
    # def eq_constraint(self):
    #     return [self.pen_ratio_cons]

    @property
    def cost_NPV(self): #net present value
        cost_NPV = self.wec.cost_NPV + self.aqua.cost_NPV + self.es.cost_NPV 
        return cost_NPV
    
    @property
    def cost_NPV_diesel(self): #net present value
        cost_NPV_diesel = self.dieselgen.cost_NPV + self.aqua.cost_NPV
        return cost_NPV_diesel
    
    @property
    def cost_per_yield(self): #net present value
       # cost_per_yield = self.cost_NPV / (self.aqua.fish_yield * self.aqua.lifetime)

        PVIF_sum = 0 # present value interest factor (PVIF)
        for i in range(self.aqua.lifetime):
            PVIF_sum += 1 / ((1+self.aqua.discount_rate)**(i+1))
        cost_per_yield = self.cost_NPV / (self.aqua.fish_yield * PVIF_sum)
        return cost_per_yield
    
    def power(self):
        self.dieselgen.power(np.max(self.aqua.power))
        return self.wec.P_electrical
    
    def plot_variable(self):
        self.aqua.fish.plot_variable()
        return
    
    def carrying_capacity_print(self):
        return self.aqua.TPF_O2, self.aqua.carrying_capacity
    
    def plot_power(self):
        fig, axes = plt.subplots(5, 1, figsize=(10, 30))

        ax1 = axes[0]
        ax1.plot(self.wec.P_mechanical)
        ax1.set(xlabel='time [hour]', ylabel='P_mechanical(kW)');
        ax1.grid()
        #plt.ylim(-10, 10)

        ax1 = axes[1]
        ax1.plot(self.wec.P_electrical)
        ax1.set(xlabel='time [hour]', ylabel='P_electrical(kW)');
        ax1.grid()
        #plt.ylim(0, 20)

        ax1 = axes[2]
        ax1.plot(self.aqua.power)
        ax1.set(xlabel='time [hour]', ylabel='pen power(kW)');
        ax1.grid()
        #plt.ylim(-10, 10)

        ax1 = axes[3]
        ax1.plot(self.es.P_diff)
        ax1.set(xlabel='time [hour]', ylabel='P_diff(kW)');
        ax1.grid()
        #plt.ylim(-20, 20)

        ax3 = axes[4]
        ax3.plot(self.es.power)
        plt.axhline(y=self.es.total_size, color='r', linestyle='-')
        plt.axhline(y=self.es.total_size * self.es.soc_uplimit, color='green', linestyle='dotted')
        plt.axhline(y=self.es.total_size * self.es.soc_downlimit, color='green', linestyle='dotted')
        ax3.set(xlabel='time [hour]', ylabel='P_stored (kW)')
        ax3.grid()
        
        plt.tight_layout()
        plt.show()


        fig, axes = plt.subplots(3, 1, figsize=(10, 18))

        ax1 = axes[0]
        ax1.plot(self.aqua.power_summer, label='power_summer total')
        ax1.plot(self.aqua.power_winter, label='power_winter total')
        ax1.set(xlabel='time [hour]', ylabel='power total')
        ax1.legend()
        ax1.grid()
        
        ax1 = axes[1]
        ax1.plot(self.aqua.summer_feedbarge_power * self.aqua.feedbarge_number, label='summer_feedbarge_power') 
        ax1.plot(self.aqua.summer_lighting_power_per_kg * self.aqua.fish_yield, label='summer_lighting_power') 
        ax1.plot(self.aqua.summer_equipment_power_per_kg * self.aqua.fish_yield, label='summer_equipment_power') 
        ax1.set(xlabel='time [hour]', ylabel='power')
        ax1.legend()
        ax1.grid()

        ax1 = axes[2]
        ax1.plot(self.aqua.winter_feedbarge_power * self.aqua.feedbarge_number, label='winter_feedbarge_power') 
        ax1.plot(self.aqua.winter_lighting_power_per_kg * self.aqua.fish_yield, label='winter_lighting_power') 
        ax1.plot(self.aqua.winter_equipment_power_per_kg * self.aqua.fish_yield, label='winter_equipment_power') 
        ax1.set(xlabel='time [hour]', ylabel='power')
        ax1.legend()
        ax1.grid()


# ============================================================================ #
#                       Set Inputs to Objects Classes                          #
# ============================================================================ #

def input_merge(x_in, x_name, p):
    # merge input dicts
    
    x_list = variable_lookup(x_name)
    x = {}
    for i in range(len(x_list)):
        x[x_list[i]] = x_in[i]
 
    valid_point = False
    gis_data = []
 
    if 'handler' in p:
        if 'pos_env' in x_name:
            gis_data =  p['handler'].query(x_in[1], x_in[0]) #import_gis_data(x_in[0], x_in[1])
        else:
            gis_data =  p['handler'].query(p['pos_long'], p['pos_lat']) #import_gis_data(x_in[0], x_in[1])
        
        duration = 8760
        p['U'] = float(gis_data["current [m/s]"])
        p['O2_in'] = float(gis_data["oxygen [mg/l]"])
        p['salinity'] = float(gis_data["salinity [PSU]"])
        #p['temp'] = float(gis_data["temperature [°C]"])
        p['bathymetry'] = (-float(gis_data["bathymetry [m]"]))
        p['distance'] = float(gis_data["distance to port [m]"]) / 1000
        #p['wave_energy_period'] = np.ones(duration) * float(gis_data["period [s]"])
        #p['wave_height'] = np.ones(duration) * float(gis_data["height [m]"])
        valid_point = gis_data["ok-conditions"].bool() and gis_data["ok-scope"].bool() #and gis_data["ok-conflicts"].bool()
    else:
        print('GIS handler is needed!')
        exit()
        

    ins = {**x, **p}
    #print(x)

    # create objects    
    wec = WEC(ins['wave_height'], ins['wave_energy_period'], 
              ins['capture_width'], ins['capture_width_ratio_dict'],
              ins['wave_damping_dict'], ins['wec_type'],
              ins['capacity_factor'], ins['eta'], ins['float_diameter'],
              ins['wec_CapEx_ref'], ins['wec_OpEx_ref'], ins['lifetime'], ins['discount_rate'], ins['P_wec_rated'])
    
    aqua = Aqua(ins['F_f'], ins['F_p'], ins['F_c'], ins['A_f'], ins['A_p'], ins['A_c'],
                ins['O_f'], ins['O_p'], ins['O_c'], ins['C_f'], ins['C_p'], ins['C_c'],
                ins['P_f'], ins['P_p'], ins['tau'], ins['loss_rate'], ins['harvest_weight'], 
                ins['O2_min'], ins['U_min'], ins['U_max'], ins['temp_min'], ins['temp_max'], 
                ins['salinity_min'], ins['salinity_max'], ins['fish_life_cycle'], ins['fingerling_weight'], ins['fingerling_unit_cost'],
                
                ins['pen_diameter'], #ins['pen_height'], 
                ins['pen_depth'], ins['stock_density'], 
                ins['num_pens'], ins['spacing'], ins['water_temp'], 
                ins['O2_in'], ins['U'], ins['salinity'], ins['permeability'], ins['bathymetry'],
                ins['pos_lat'], ins['pos_long'], 
                ins['pen_netting_CapEx_ref'], ins['pen_struct_CapEx_ref'], ins['feedbarge_CapEx_ref'], ins['feedbarge_OpEx_ref'], 
                ins['lifetime'], ins['discount_rate'],
                ins['FCR'], ins['feed_unit_cost'], ins['feedbarge_unit_capacity'], ins['feedbarge_unit_feedlines'],
                ins['summer_feedbarge_power'], ins['summer_lighting_power_per_kg'], ins['summer_equipment_power_per_kg'],
                ins['winter_feedbarge_power'], ins['winter_lighting_power_per_kg'], ins['winter_equipment_power_per_kg'],
                
                ins['vessel_fuel_consump_rate'], ins['vessel_fuel_cost'], 
                ins['captain_salary'], ins['crew_salary'], ins['crew_num'], 
                ins['time_feed'], ins['vessel_velocity'], ins['travel_number'], ins['distance']
                )
    
    es = ES(ins['es_eta'], ins['es_CapEx_ref'], ins['es_OpEx_ref'], ins['lifetime'], ins['discount_rate'], ins['es_soc_uplimit'], ins['es_soc_downlimit'])

    dieselgen = DieselGen(ins['diesel_fuel_consump_rate'], ins['diesel_fuel_cost'], ins['diesel_eta'], ins['diesel_load_level'],
                          ins['diesel_CapEx_ref'], ins['diesel_OpEx_ref'], ins['lifetime'], ins['discount_rate'])

    return wec, aqua, es, dieselgen, valid_point, gis_data

# ============================================================================ #
#                       Define DVs and Parameters                              #
# ============================================================================ #

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
        # var_list.append('pen_height')
        var_list.append('stock_density')

    if any('x_disc_pen' in i for i in var_category_names):
        var_list.append('num_pens')

    if any('p_pen' in i for i in var_category_names):
        var_list.append('spacing')
        var_list.append('pen_depth')
        var_list.append('permeability')
        var_list.append('pen_netting_CapEx_ref')
        var_list.append('pen_struct_CapEx_ref')
    
    if any('p_feedbarge' in i for i in var_category_names):
        var_list.append('feedbarge_CapEx_ref')
        var_list.append('feedbarge_OpEx_ref')
        var_list.append('feedbarge_unit_capacity')
        var_list.append('feedbarge_unit_feedlines')
    
    if any('p_pen_power' in i for i in var_category_names):
        var_list.append('summer_feedbarge_power')
        var_list.append('summer_lighting_power_per_kg')
        var_list.append('summer_equipment_power_per_kg')
        var_list.append('winter_feedbarge_power')
        var_list.append('winter_lighting_power_per_kg')
        var_list.append('winter_equipment_power_per_kg')

    if any('pos_env' in i for i in var_category_names):
        var_list.append('pos_lat')
        var_list.append('pos_long')
    
    if any('gis_handler' in i for i in var_category_names):
        var_list.append('handler')
        
    if any('x_env' in i for i in var_category_names):
        #var_list.append('temp')
        var_list.append('O2_in')
    
    if any('x_wave_env' in i for i in var_category_names):
        var_list.append('wave_height')
        var_list.append('wave_energy_period')
        var_list.append('water_temp')
        
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
        var_list.append('wec_CapEx_ref')
        var_list.append('wec_OpEx_ref')
        var_list.append('lifetime')
        var_list.append('discount_rate')
        var_list.append('P_wec_rated')
    
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
        var_list.append('fish_life_cycle') 
        var_list.append('fingerling_weight') 
        var_list.append('fingerling_unit_cost') 
    
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

    if any('p_es' in i for i in var_category_names):
        var_list.append('es_eta')
        var_list.append('es_CapEx_ref')
        var_list.append('es_OpEx_ref')
        var_list.append('es_soc_uplimit')
        var_list.append('es_soc_downlimit')

    if any('p_diesel' in i for i in var_category_names):
        var_list.append('diesel_eta')
        var_list.append('diesel_fuel_consump_rate')
        var_list.append('diesel_fuel_cost')
        var_list.append('diesel_load_level')
        var_list.append('diesel_CapEx_ref')
        var_list.append('diesel_OpEx_ref')

    if len(var_list)==0:
        print('Your input did not match any of the category names.', var_category_names)
    

    return var_list

# ============================================================================ #
#            Define Values for Parameters and Initial Values for DVs           #
# ============================================================================ #

def default_values(var_category_names):
    vals = {}
    wec_types = (['attenuator','terminator','point absorber'], '[-]')
    capture_width_ratios = ([0.16, 0.34, 0.16], '[-]') 
    wave_dampings = ([0, 0.13, 0.17], '[-]')            

    if any('x_wec' in i for i in var_category_names):
        vals['capture_width'] = (20, '[m]') #20

    if any('x_type_wec' in i for i in var_category_names):
        vals['wec_type'] = ('point absorber', '[-]')
        
    if any('x_pen' in i for i in var_category_names):
        vals['pen_diameter'] = (30, '[m]') 
        # vals['pen_height'] = (15, '[m]')  
        vals['stock_density'] = (20 , '[kg/m^3]')
   
    if any('x_disc_pen' in i for i in var_category_names):
        vals['num_pens'] = (12, '[-]')  #{5, 12, 40}

    if any('p_pen' in i for i in var_category_names):
        vals['spacing'] = (150, '[m]')
        vals['pen_depth'] = (10, '[m]')   
        vals['permeability'] = (0.8, '[-]')   
        #vals['pen_CapEx_ref'] = (100, '[$/m^3]')    # 80 $/m^3 for net pen + 20 $/m^3 for mooring
        vals['pen_netting_CapEx_ref'] = (59, '[$/m^2]')
        vals['pen_struct_CapEx_ref'] = (817, '[$/m]')
    
    if any('p_feedbarge' in i for i in var_category_names):
        vals['feedbarge_CapEx_ref'] = (1867090, '[$]')
        vals['feedbarge_OpEx_ref'] = (0, '[$]') 
        vals['feedbarge_unit_capacity'] =  (200000, '[kg]') 
        vals['feedbarge_unit_feedlines'] = (6, '[-]') 
    
    if any('p_pen_power' in i for i in var_category_names):
        vals['summer_feedbarge_power'] = ({}, '[kw]')
        vals['summer_lighting_power_per_kg'] = ({}, '[kw/kg]')
        vals['summer_equipment_power_per_kg'] = ({}, '[kw/kg]')
        vals['winter_feedbarge_power'] = ({}, '[kw]')
        vals['winter_lighting_power_per_kg'] = ({}, '[kw/kg]')
        vals['winter_equipment_power_per_kg'] = ({}, '[kw/kg]')

    if any('pos_env' in i for i in var_category_names):
        vals['pos_lat'] = (44.103, 'm') #42.0
        vals['pos_long'] = (-68.112, 'm') #-70.0
    
    if any('gis_handler' in i for i in var_category_names):
        vals['handler'] = ({}, '[-]') 

    if any('x_env' in i for i in var_category_names):
        #vals['temp'] = (9.53, 'C') #10.29
        vals['O2_in'] = (9.8289,'[mg/l]') #9.5
    
    if any('x_wave_env' in i for i in var_category_names):
        vals['wave_height'] = (1.37, '[m]')  
        vals['wave_energy_period'] = (6.44, '[s]')
        vals['water_temp'] = (9.53, 'C') #10.29
    
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
        vals['wec_CapEx_ref'] = (5142789, '[$]') #cost in 2023
        vals['wec_OpEx_ref'] = (123477, '[$]')   #cost in 2023
        vals['lifetime'] = (20, '[year]')
        vals['discount_rate'] = (0.07, '[-]')
        vals['P_wec_rated'] = (300, '[kW]')
        
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
        vals['harvest_weight'] = (4.6, '[kg/fish]')     #Fish Harvest Size [kg/fish]  # can be between 4 [kg] to 6 [kg] Value is now coming from fish growth model
        vals['O2_min'] = (4.41, '[mg/l]')               #Dissolved Oxygen Threshold [%]
        vals['U_min'] = (0.01,'[m/s]')                  # From NorthEast U.S. Environment
        vals['U_max'] = (2,'[m/s]')
        vals['temp_min'] = (2,'[C]')
        vals['temp_max'] = (20,'[C]')
        vals['salinity_min'] = (30,'[PSU]')
        vals['salinity_max'] = (35,'[PSU]')
        vals['FCR'] = (1.35,'[kgFeed/kgFish]')          #Feed Conversion Ratio [kgFeed/kgFish]
        vals['feed_unit_cost'] = (1.54,'[$/kgFeed]')   #[$/kgFeed]
        vals['fish_life_cycle'] = (365,'[day]')        #[day] 
        vals['fingerling_weight'] = (.200,'[kg]')      # [kg] previously 0.1 kg
        vals['fingerling_unit_cost'] = (2.5,'[$/smolt]')    # cost per each smolt between 85 gr to 150 gr.
    
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
    
    if any('p_es' in i for i in var_category_names):
        vals['es_eta'] = (0.92,'[-]')
        vals['es_soc_uplimit'] = (.85,'%')  #0.85
        vals['es_soc_downlimit'] = (.15,'%') #.15
        vals['es_CapEx_ref'] = (271 * 1.14, '[$/kWh]') #change to $ in 2023   * 1.14
        vals['es_OpEx_ref'] = (0, '[$/kWh]')

    if any('p_diesel' in i for i in var_category_names):
        vals['diesel_eta'] = (.40,'[%]')
        vals['diesel_fuel_consump_rate'] = (7.6,'[gal/h]')  # for 135 kW
        vals['diesel_fuel_cost'] = (3.5,'[$/gal]')
        vals['diesel_load_level'] = (0.75,'[%]') 
        vals['diesel_CapEx_ref'] = (0.3*1.14,'[$/kWh]')
        vals['diesel_OpEx_ref'] = (0,'[$/kWh]')

    #assert(fieldnames(vals) == variable_lookup(var_category_names));
    return vals

# ============================================================================ #
#                       Define Bounds for DVs                                  #
# ============================================================================ #

def bnds_values(var_category_names):
    bnds = {}

    if any('x_wec' in i for i in var_category_names):
        bnds['capture_width'] = (1, 100)      #[m] (1, 40)
    
    if any('x_pen' in i for i in var_category_names):
        bnds['pen_diameter'] = (10, 45)      #[m]  (10, 45)
        # bnds['pen_height'] = (10, 45)        #[m]   (10, 30)
        bnds['stock_density'] = (1, 20)     #[kg/m^3]  (10, 20)
    
    if any('x_disc_pen' in i for i in var_category_names):
        bnds['num_pens'] = (5, 40)         #[-] (5, 40)

    if any('pos_env' in i for i in var_category_names):
        bnds['pos_lat'] = (38.4, 45.2)        #[m]
        bnds['pos_long'] = (-75.8, -65.7)     #[m]
    
    if any('x_env' in i for i in var_category_names):
        #bnds['temp'] = (1, 50)              #[C] 
        bnds['O2_in'] = (1, 50)             #[mg/l]
        bnds['wave_height'] = (0.2, 3)      #[m]
        bnds['wave_period'] = (1, 12)       #[s]
    
    if any('x_es' in i for i in var_category_names):
        bnds['es_size'] = (5, 5000)          #[kWh] (5,500) without cost travel from shore

    return bnds
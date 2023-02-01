import numpy as np
from objects import *    

def obj(x_in, x_name, p_in: dict):
    cost_per_yield, _, _, _, _ = obj_terms(x_in, x_name, p_in)
    J = np.array(cost_per_yield)
    return J

def obj_terms(x_in, x_name, p_in: dict):
    # merge input dicts
    wec, wave_in, pen, fish = input_merge(x_in, x_name, p_in)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish_yield_func(wave_out, pen, fish)
    wec.P_gen = power(wec, wave_in)
    
    price = econ(wec, pen)
   
    # outputs
    cost_per_yield = price/fish_yield 
    return cost_per_yield, price, fish_yield, pen.price, wec.price

def ineq_constraint(x_in, x_name, p):
    # merge input dicts
    wec, wave_in, pen, fish = input_merge(x_in, x_name, p)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish_yield_func(wave_out, pen, fish)
    wec.P_gen = power(wec, wave_in)
    carrying_capacity = carrying_capacity_func(pen, fish)
    
    # power supply constraint to ensure supply power demand of net pen
    P_gen_cons = (wec.annual_energy - pen.power) / wec.annual_energy

    # fish yield constraint to ensure a healthy offshore environment
    fish_yield_cons = (carrying_capacity - fish_yield) / carrying_capacity
 
    # net pen geometry constraint to present a practical design and ratio between diameter and height
    pen_ratio_low_cons = (pen.D - pen.H) / pen.D
    pen_ratio_up_cons = (3*pen.H - pen.D) / (3*pen.H)
    
    # outputs
    g = np.array([P_gen_cons, fish_yield_cons, pen_ratio_low_cons, pen_ratio_up_cons])

    return g


def eq_constraint(x_in, x_name, p):
    h = np.array([0])
        
    return h

def input_merge(x_in, x_name, p):
    # merge input dicts
    
    x_list = variable_lookup(x_name)
    x = {}
    for i in range(len(x_list)):
        x[x_list[i]] = x_in[i]
 
    ins = {**x, **p}

    # create objects
    wave_in = Wave(ins['wave_height'], ins['wave_period'])
    
    wec = WEC(ins['capture_width'], ins['capture_width_ratio_dict'],
            ins['wave_damping_dict'], ins['wec_type'], ins['wec_unit_cost'])

    pen = Pen(ins['pen_diameter'], ins['pen_height'], ins['pen_depth'], ins['stock_density'], 
              ins['num_pens'], ins['spacing'], ins['pen_unit_cost'], ins['temp'], 
              ins['O2_in'], ins['U'], ins['salinity'], ins['permeability'])
    
    fish = Fish(ins['F_f'], ins['F_p'], ins['F_c'], ins['A_f'], ins['A_p'], ins['A_c'],
                ins['O_f'], ins['O_p'], ins['O_c'], ins['C_f'], ins['C_p'], ins['C_c'],
                ins['P_f'], ins['P_p'], ins['tau'], ins['loss_rate'], ins['harvest_weight'], 
                ins['O2_min'], ins['U_min'], ins['U_max'], ins['temp_min'], ins['temp_max'], 
                ins['salinity_min'], ins['salinity_max'])
                
    return wec, wave_in, pen, fish

def plot_variable(x_in, x_name, p):
    # merge input dicts
    wec, wave_in, pen, fish = input_merge(x_in, x_name, p)
    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish_yield_func(wave_out,pen, fish)
    wec.P_gen = power(wec, wave_in)
    carrying_capacity = carrying_capacity_func(pen, fish)
    P_gen_cons = wec.annual_energy - pen.power
    fish_yield_cons = carrying_capacity - fish_yield
    fish.plot_variable
    return

def power(wec: WEC, wave: Wave) -> float:
    assert(isinstance(wec,WEC))
    assert(isinstance(wave,Wave))
    P_gen = wave.power * wec.capture_width * wec.capture_width_ratio
    
    return P_gen


def P_rated(x_in, x_name, p_in: dict):
    # merge input dicts
    wec, wave_in, pen, fish = input_merge(x_in, x_name, p_in)
    wec.P_gen = power(wec, wave_in)
    
    return wec.P_gen/wec.capacity_factor
    

def carrying_capacity_print(x_in, x_name, p_in: dict):
    # merge input dicts
    wec, wave_in, pen, fish = input_merge(x_in, x_name, p_in)
    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish_yield_func(wave_out, pen, fish)
    wec.P_gen = power(wec, wave_in)
    carrying_capacity = carrying_capacity_func(pen, fish)
    
    return pen.TPF_O2, carrying_capacity

def econ(wec: WEC, pen: Pen) -> float:
    assert(isinstance(wec,WEC))
    assert(isinstance(pen,Pen))
    price = wec.price + pen.price
    return price

def wave_climate(wec: WEC, wave: Wave) -> Wave:
    assert(isinstance(wec,WEC))
    assert(isinstance(wave,Wave))
    Hs = wec.wave_damping * wave.Hs
    T = wave.T
    wave_out = Wave(Hs,T)
    return wave_out

def fish_yield_func(wave: Wave, pen: Pen, fish: Fish) -> float:
    assert(isinstance(wave,Wave))
    assert(isinstance(pen,Pen))
    survival_rate = (1-fish.loss_rate)
    pen.fish_yield = pen.n * pen.SD * survival_rate * pen.volume
    return pen.fish_yield

def carrying_capacity_func(pen: Pen, fish: Fish) -> float:
    assert(isinstance(pen,Pen))
    assert(isinstance(fish,Fish))
    return pen.carrying_capacity(fish)


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
        
    if any('x_env' in i for i in var_category_names):
        var_list.append('temp')
        var_list.append('O2_in')
        var_list.append('wave_height')
        var_list.append('wave_period')
        
    if any('p_env' in i for i in var_category_names):
        var_list.append('salinity')
        var_list.append('U')
        
    if any('p_wec' in i for i in var_category_names):
        var_list.append('wec_unit_cost')
        var_list.append('capture_width_ratio_dict')
        var_list.append('wave_damping_dict')
    
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

    if len(var_list)==0:
        print('Your input did not match any of the category names.', var_category_names)
    

    return var_list


def default_values(var_category_names):
    vals = {}
    wec_types = (['attenuator','terminator','point absorber'], '[-]')
    capture_width_ratios = ([0.16, 0.34, 0.35], '[-]') 
    wave_dampings = ([0, 0.13, 0.17], '[-]')            

    if any('x_wec' in i for i in var_category_names):
        vals['capture_width'] = (76, '[m]')

    if any('x_type_wec' in i for i in var_category_names):
        vals['wec_type'] = ('point absorber', '[-]')
        
    if any('x_pen' in i for i in var_category_names):
        vals['pen_diameter'] = (25, '[m]') 
        vals['pen_height'] = (10, '[m]')  
        vals['stock_density'] = (30 , '[kg/m^3]')
   
    if any('p_pen' in i for i in var_category_names):
        vals['num_pens'] = (18, '[-]')  
        vals['spacing'] = (150, '[m]') 
        vals['pen_depth'] = (10, '[m]')   
        vals['pen_unit_cost'] = (100, '[$/m^3]')    # 80 $/m^3 for net pen + 20 $/m^3 for mooring
        vals['permeability'] = (0.8, '[-]')      
        
    if any('x_env' in i for i in var_category_names):
        vals['temp'] = (10.07, 'C')
        vals['O2_in'] = (9.44,'[mg/l]')
        vals['wave_height'] = (1.04, '[m]')
        vals['wave_period'] = (5.73, '[s]') 
    
    if any('p_env' in i for i in var_category_names):
        vals['salinity'] = (31.74, '[PSU]')
        vals['U'] = (.1, '[m/s]')
    
    if any('p_wec' in i for i in var_category_names):
        vals['wec_unit_cost'] = (0.45 * 1.19, '[$/kWh]')   # 'point absorber' * inflation rate from 2014 to 2022
        vals['capture_width_ratio_dict'] = (dict(zip(wec_types[0], capture_width_ratios[0])), '[-]')
        vals['wave_damping_dict'] = (dict(zip(wec_types[0], wave_dampings[0])), '[-]')
        
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
        vals['temp_max'] = (28,'[C]')
        vals['salinity_min'] = (30,'[PSU]')
        vals['salinity_max'] = (35,'[PSU]')
    
    #assert(fieldnames(vals) == variable_lookup(var_category_names));
    return vals

def bnds_values(var_category_names):
    bnds = {}

    if any('x_wec' in i for i in var_category_names):
        bnds['capture_width'] = (20, 85)     #[m]
    
    if any('x_pen' in i for i in var_category_names):
        bnds['pen_diameter'] = (25, 45)       #[m] 
        bnds['pen_height'] = (10, 30)        #[m]
        bnds['stock_density'] = (15, 30)     #[kg/m^3]
    
    if any('x_env' in i for i in var_category_names):
        bnds['temp'] = (7, 22)              #[C]
        bnds['O2_in'] = (7, 10)             #[mg/l]
        bnds['wave_height'] = (0.2, 3)      #[m]
        bnds['wave_period'] = (1, 12)       #[s]
    
    return bnds
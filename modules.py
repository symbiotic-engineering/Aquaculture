import numpy as np
from objects import *
from typing import Tuple
    

def obj(x_in, x_name, p_in: dict):
    cost_per_yield, price, fish_yield = obj_terms(x_in, x_name, p_in)
    J = np.array(cost_per_yield)
    return J

def obj_terms(x_in, x_name, p_in: dict):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, x_name, p_in)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish(wave_out,pen)
    wec.P_gen = power(wec, wave_in)
    
    price = econ(wec, pen)
   
    # outputs
    cost_per_yield = price/fish_yield 
    #print(x_in, price, fish_yield, cost_per_yield)
    
    return cost_per_yield, price, fish_yield

def ineq_constraint(x_in, x_name, p):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, x_name, p)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish(wave_out,pen)
    wec.P_gen = power(wec, wave_in)
    carrying_capacity = environment(pen)
    
    P_gen_cons = wec.annual_energy - pen.power
    #P_gen_cons = wec.P_gen - pen.power
    fish_yield_cons = carrying_capacity - fish_yield
    #print(x_in, carrying_capacity, fish_yield, wec.annual_energy, pen.power)

    # outputs
    g = np.array([P_gen_cons, fish_yield_cons])

    return g


def eq_constraint(x_in, x_name, p):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, x_name, p)
    
    h = np.array([])
        
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
            ins['num_pens'], ins['spacing'], ins['pen_unit_cost'], ins['loss_rate'],
            ins['harvest_weight'], ins['temp'], 
            ins['O2_in'],ins['O2_min'],ins['P_f'],ins['P_p'],ins['U_min'],
            ins['tau'],ins['permeability'],ins['F_f'],ins['F_p'],
            ins['F_c'],ins['A_f'],ins['A_p'],ins['A_c'],ins['O_f'],ins['O_p'],
            ins['O_c'],ins['C_f'],ins['C_p'],ins['C_c'])
    
    return wec, wave_in, pen

def plot_variable(x_in, x_name, p):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, x_name, p)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish(wave_out,pen)
    wec.P_gen = power(wec, wave_in)
    carrying_capacity = environment(pen)
    
    P_gen_cons = wec.annual_energy - pen.power
    fish_yield_cons = carrying_capacity - fish_yield
    pen.plot_variable
    return

def power(wec: WEC, wave: Wave) -> float:
    assert(isinstance(wec,WEC))
    assert(isinstance(wave,Wave))
    
    P_gen = wave.power * wec.capture_width * wec.capture_width_ratio
    
    return P_gen

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

def fish(wave: Wave, pen: Pen) -> float:
    assert(isinstance(wave,Wave))
    assert(isinstance(pen,Pen))
              
    if wave.Hs > 3:
        extra_loss_rate = 0.1
    else:
        extra_loss_rate = 0

    survival_rate = (1-(pen.loss_rate + extra_loss_rate))
    pen.fish_yield = pen.n * pen.SD * survival_rate * pen.volume * pen.harvest_weight
    return pen.fish_yield

def environment(pen: Pen) -> float:
    assert(isinstance(pen,Pen))

    return pen.carrying_capacity


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

    if any('p_pen' in i for i in var_category_names):
        var_list.append('num_pens')
        var_list.append('pen_unit_cost')
        var_list.append('permeability')
        
    if any('x_env' in i for i in var_category_names):
        var_list.append('temp')
        var_list.append('O2_in')
        var_list.append('salinity')
        var_list.append('wave_height')
        var_list.append('wave_period')
        var_list.append('U_min')
        
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
        var_list.append('O2_min')          #Dissolved Oxygen Threshold [%]
        var_list.append('tau')             #Inverse Temperature Scale [1/C]
        var_list.append('loss_rate')       #Fish Loss Rate [-]
        var_list.append('harvest_weight')  #Fish Harvest Size [kg/fish]
    
    if len(var_list)==0:
        print('Your input did not match any of the category names.', var_category_names)
    

    return var_list


def default_values(var_category_names):
    vals = {}
    wec_types = (['attenuator','terminator','point absorber'], '[-]')
    capture_width_ratios = ([0.16, 0.34, 0.35], '[-]') 
    wave_dampings = ([0, 0.13, 0.17], '[-]')            

    if any('x_wec' in i for i in var_category_names):
        vals['capture_width'] = (30, '[m]')     

    if any('x_type_wec' in i for i in var_category_names):
        vals['wec_type'] = ('point absorber', '[-]')
        
    if any('x_pen' in i for i in var_category_names):
        vals['pen_diameter'] = (30, '[m]')      
        vals['pen_height'] = (15, '[m]')        
        vals['spacing'] = (150, '[m]')          
        vals['stock_density'] = (10 , '[kg/m^3]') 
        vals['pen_depth'] = (80, '[m]')         
   
    if any('p_pen' in i for i in var_category_names):
        vals['num_pens'] = (18, '[-]')          
        
    if any('x_env' in i for i in var_category_names):
        vals['temp'] = (16, 'C')               
        vals['salinity'] = (33, '[PSU]')        
        vals['O2_in'] = (8,'[mg/l]')            
        vals['U_min'] = (0.25,'[m/s]')          
        vals['wave_height'] = (2.65, '[m]')     
        vals['wave_period'] = (8.33, '[s]')     
    
    if any('p_wec' in i for i in var_category_names):
        vals['wec_unit_cost'] = (0.45, '[$/kWh]')   # 'point absorber'
        vals['pen_unit_cost'] = (100, '[$/m^3]')    # 80 $/m^3 for net pen + 20 $/m^3 for mooring
        vals['permeability'] = (0.8, '[-]')
        vals['capture_width_ratio_dict'] = (dict(zip(wec_types[0], capture_width_ratios[0])), '[-]')
        vals['wave_damping_dict'] = (dict(zip(wec_types[0], wave_dampings[0])), '[-]')
        
    if any('p_fish_salmon' in i for i in var_category_names):
        vals['F_p'] = (0.45, '[-]')
        vals['F_c'] = (0.07, '[-]')
        vals['F_f'] = (0.15, '[-]')
        vals['A_p'] = (0.97, '[-]')
        vals['A_c'] = (0.6, '[-]')
        vals['A_f'] = (0.9, '[-]')
        vals['O_p'] = (1.89, '[gO2/gProtein]')
        vals['O_c'] = (1.07, '[gO2/gCarbohydrate]')
        vals['O_f'] = (2.91, '[gO2/gFat]')
        vals['C_p'] = (5650, '[cal/g]')
        vals['C_c'] = (4100, '[cal/g]')
        vals['C_f'] = (9450, '[cal/g]')
        vals['P_f'] = (0.18, '[-]')
        vals['P_p'] = (0.18, '[-]')
        vals['O2_min'] = (0.9, '[%]')
        vals['tau'] = (0.08, '[1/C]')
        vals['loss_rate'] = (0.15, '[-]')
        vals['harvest_weight'] = (4, '[kg/fish]')
    
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
        bnds['capture_width'] = (1, 40)     #[m]
    
    if any('x_pen' in i for i in var_category_names):
        bnds['pen_diameter'] = (3, 50)      #[m]
        bnds['pen_height'] = (5, 30)        #[m]
        bnds['spacing'] = (50, 300)         #[m]
        bnds['stock_density'] = (1, 30)     #[kg/m^3] 50
        bnds['pen_depth'] = (40, 120)       #[m] (40, 120)
    
    if any('x_env' in i for i in var_category_names):
        bnds['temp'] = (7, 22)              #[C]
        bnds['salinity'] = (26, 36)         #[PSU]
        bnds['O2_in'] = (7, 10)             #[mg/l]
        bnds['U_min'] = (0.01, 0.55)        #[m/s]
        bnds['wave_height'] = (1, 5)        #[m]
        bnds['wave_period'] = (1, 12)       #[s]
    
    return bnds
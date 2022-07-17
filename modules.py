import numpy as np
from objects import *
from typing import Tuple

def obj(x_in, x_name, p_in: dict):
    # merge input dicts
    #print(x_in, x_name, p_in)
    wec, wave_in, pen = input_merge(x_in, x_name, p_in)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish(wave_out,pen)
    price = econ(wec, pen)

    # outputs
    cost_per_yield = price/fish_yield 
    J = np.array(cost_per_yield)
    
    return J


def ineq_constraint(x_in, x_name, p):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, x_name, p)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    pow = power(wec, wave_in)
    carrying_capacity = environment(pen)

    # outputs
    g = np.array([pow, carrying_capacity])

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
    wec = WEC(ins['capture_width'], ins['capture_width_ratio_dict'],
            ins['wave_damping_dict'], ins['wec_type'], ins['wec_unit_cost'])

    wave_in = Wave(ins['wave_height'], ins['wave_period'])

    pen = Pen(ins['pen_diameter'], ins['pen_height'], ins['stocking_density'], 
            ins['num_pens'], ins['spacing'], ins['pen_unit_cost'], ins['loss_rate'],
            ins['harvest_weight'], ins['temp'], 
            ins['O2_in'],ins['O2_min'],ins['P_f'],ins['P_p'],ins['U_min'],
            ins['tau'],ins['permeability'],ins['F_f'],ins['F_p'],
            ins['F_c'],ins['A_f'],ins['A_p'],ins['A_c'],ins['O_f'],ins['O_p'],
            ins['O_c'],ins['C_f'],ins['C_p'],ins['C_c'])
    
    return wec, wave_in, pen

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
    fish_yield = pen.SD * survival_rate * pen.volume * pen.harvest_weight

    return fish_yield

def environment(pen: Pen) -> float:
    assert(isinstance(pen,Pen))

    return pen.carrying_capacity


def variable_lookup(var_category_names):
    # input is a cell array containing some subset of the following categories:
    # x_wec, x_wec_type, x_pen, p_pen, p_env, p_consts, p_fish
    # output is a cell array with variable names matching the input categoeies.

    var_list = []
    if any('x_wec' in i for i in var_category_names):
        var_list.append('capture_width')

    if any('x_wec_type' in i for i in var_category_names):
        var_list.append('wec_type')
        
    if any('x_pen' in i for i in var_category_names):
        var_list.append('pen_diameter')
        var_list.append('pen_height')

    if any('p_pen' in i for i in var_category_names):
        var_list.append('num_pens')
        var_list.append('spacing')
        var_list.append('stocking_density')
        
    if any('x_env' in i for i in var_category_names):
        var_list.append('temp')
        var_list.append('O2_in')
        var_list.append('min_current_speed')
        var_list.append('wave_height')
        var_list.append('wave_period')
        var_list.append('U_min')
    
    if any('p_consts' in i for i in var_category_names):
        var_list.append('wec_unit_cost')
        var_list.append('pen_unit_cost')
        var_list.append('permeability')
        var_list.append('capture_width_ratio_dict')
        var_list.append('wave_damping_dict')
    
    if any('p_fish' in i for i in var_category_names):
        var_list.append('F_f')
        var_list.append('F_p')
        var_list.append('F_c')
        var_list.append('A_f')
        var_list.append('A_p')
        var_list.append('A_c')
        var_list.append('O_f')
        var_list.append('O_p')
        var_list.append('O_c')
        var_list.append('C_f')
        var_list.append('C_p')
        var_list.append('C_c')
        var_list.append('P_f')
        var_list.append('P_p')
        var_list.append('O2_min')
        var_list.append('tau')
        var_list.append('loss_rate')
        var_list.append('harvest_weight')
    
    if len(var_list)==0:
        print('Your input did not match any of the category names.', var_category_names)
    

    return var_list


def default_values(var_category_names):
    vals = {}
    wec_types = ['attenuator','terminator','point absorber']
    capture_width_ratios = [0.5, 0.5, 0.5]
    wave_dampings = [0.5, 0.5, 0.5]

    if any('x_wec' in i for i in var_category_names):
        vals['capture_width'] = 10

    if any('x_wec_type' in i for i in var_category_names):
        vals['wec_type'] = 'point absorber'
        
    if any('x_pen' in i for i in var_category_names):
        vals['pen_diameter'] = 1
        vals['pen_height'] = 0.5
   
    if any('p_pen' in i for i in var_category_names):
        vals['num_pens'] = 8
        vals['spacing'] = 10
        vals['stocking_density'] = 1
        
    if any('x_env' in i for i in var_category_names):
        vals['temp'] = 1
        vals['salinity'] = 1
        vals['O2_in'] = 1
        vals['min_current_speed'] = 1
        vals['wave_height'] = 2
        vals['wave_period'] = 5
        vals['U_min'] = 1
    
    if any('p_consts' in i for i in var_category_names):
        vals['wec_unit_cost'] = 1000
        vals['pen_unit_cost'] = 1000
        vals['permeability'] = 1
        vals['capture_width_ratio_dict'] = dict(zip(wec_types, capture_width_ratios))
        vals['wave_damping_dict'] = dict(zip(wec_types, wave_dampings))
    
    if any('p_fish' in i for i in var_category_names):
        vals['F_p'] = 1
        vals['F_c'] = 1
        vals['F_f'] = 1
        vals['A_p'] = 1
        vals['A_c'] = 1
        vals['A_f'] = 1
        vals['O_p'] = 1
        vals['O_c'] = 1
        vals['O_f'] = 1
        vals['C_p'] = 1
        vals['C_c'] = 1
        vals['C_f'] = 1
        vals['P_f'] = 1
        vals['P_p'] = 1
        vals['O2_min'] = 0.9
        vals['tau'] = 1
        vals['loss_rate'] = 0.15
        vals['harvest_weight'] = 1

    #assert(fieldnames(vals) == variable_lookup(var_category_names));
    return vals
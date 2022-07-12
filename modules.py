import numpy as np
from objects import *
from typing import Tuple

def obj(x_in: dict, p: dict):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, p)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish(wave_out,pen)
    price = econ(wec, pen)

    # outputs
    cost_per_yield = price/fish_yield 
    J = np.array(cost_per_yield)
    
    return J


def ineq_constraint(x_in, p):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, p)

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    pow = power(wec, wave_in)
    carrying_capacity = environment(pen)

    # outputs
    g = np.array([pow, carrying_capacity])

    return g


def eq_constraint(x_in, p):
    # merge input dicts
    wec, wave_in, pen = input_merge(x_in, p)
    
    h = np.array([])
        
    return h

def input_merge(x_in: dict, p: dict):
    # merge input dicts
    
    x = {'capture_width': x_in[0],
    'pen_diameter': x_in[1],
    'pen_height': x_in[2]}
    
    ins = {**x, **p}

    # create objects
    wec = WEC(ins['capture_width'], ins['capture_width_ratio_dict'],
            ins['wave_damping_dict'], ins['wec_type'], ins['wec_unit_cost'])

    wave_in = Wave(ins['wave_height'], ins['wave_period'])

    pen = Pen(ins['pen_diameter'], ins['pen_height'], ins['stocking_density'], 
            ins['num_pens'], ins['spacing'], ins['pen_unit_cost'], ins['loss_rate'],
            ins['harvest_weight'], ins['env_params'])
    
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
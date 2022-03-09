import numpy as np
from objects import *
from typing import Tuple

def sim(x: dict, p: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # merge input dicts
    ins = x | p

    # create objects
    wec = WEC(ins['capture_width'], ins['capture_width_ratio_dict'],
            ins['wave_damping_dict'], ins['wec_type'], ins['wec_unit_cost'])
    wave_in = Wave(ins['wave_height'], ins['wave_period'])
    pen = Pen(ins['pen_diameter'], ins['pen_height'], ins['pen_unit_cost'])

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish(wave_out,pen)
    pow = power(wec, wave_in)
    carrying_capacity = environment(pen)
    price = econ(wec, pen)

    # outputs
    cost_per_yield = price/fish_yield 
    J = np.array(cost_per_yield)
    g = np.array([pow, carrying_capacity])
    h = np.array([])

    return J, g, h


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

    if wave.Hs < 1:
        fish_yield = 1
    else:
        fish_yield = 0.5
    return fish_yield

def environment(pen: Pen) -> float:
    assert(isinstance(pen,Pen))
    
    carrying_capacity = 1
    return carrying_capacity
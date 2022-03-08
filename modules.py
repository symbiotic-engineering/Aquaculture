from objects import *

def sim(x,p):
    # merge input dicts
    ins = x | p
    
    # create objects
    wec = WEC(ins['capture_width'], ins['capture_width_ratio'], ins['wec_type'])
    wave_in = Wave(ins['wave_height'], ins['wave_period'])
    pen = Pen(ins['pen_diameter'], ins['pen_height'])

    # run each module 
    wave_out = wave_climate(wec,wave_in)
    fish_yield = fish(wave_out,pen)
    pow = power(wec, wave_in)
    carrying_capacity = environment(pen)
    price = econ(wec, pen)

    # outputs
    cost_per_yield = price/fish_yield 
    return pow, cost_per_yield, carrying_capacity


def power(wec:WEC, Wave:Wave):
    P_gen = Wave.power * wec.capture_width * wec.capture_width_ratio
    return P_gen

def econ(wec:WEC, pen:Pen):
    price = wec.price + pen.price
    return price

def wave_climate(wec:WEC, wave:Wave):
    Hs = wec.wave_damping * wave.Hs
    T = wave.T
    wave_out = Wave(Hs,T)
    return wave_out

def fish(wave:Wave, pen:Pen):
    if wave.Hs < 1:
        fish_yield = 1
    else:
        fish_yield = 0.5
    return fish_yield

def environment(pen:Pen):
    carrying_capacity = 1
    return carrying_capacity
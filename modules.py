from objects import *

def sim(wec:WEC, waveIn:Wave, pen:Pen):
    waveOut = wave_climate(wec,waveIn)
    fish_yield = fish(waveOut,pen)
    pow = power(wec, waveIn)
    carrying_capacity = environment(pen)
    price = econ(wec, pen)
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
    waveOut = Wave(Hs,T)
    return waveOut

def fish(wave:Wave, pen:Pen):
    if wave.Hs < 1:
        fish_yield = 1
    else:
        fish_yield = 0.5
    return fish_yield

def environment(pen:Pen):
    carrying_capacity = 1
    return carrying_capacity
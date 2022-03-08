from math import pi

class WEC:
    def __init__(self, capture_width, capture_width_ratio, type):
        self.capture_width = capture_width
        self.capture_width_ratio = capture_width_ratio
        self.type = type

    @property
    def price(self):
        price = self.capture_width * 1000
        return price

    @property
    def wave_damping(self):
        if self.type == 'attenuator':
            damping = .5
        elif self.type == 'point absorber':
            damping = 1
        elif self.type == 'terminator':
            damping = .75
        else: 
            print('Error: Wrong WEC type')
        return damping

class Wave:
    def __init__(self,Hs,T):
        self.Hs = Hs
        self.T = T
        self.rho = 1000
        self.g = 9.81
    
    @property
    def power(self):
        P_wave = 1/32 * 1/pi * self.rho * self.g**2 * self.Hs**2 * self.T
        return P_wave

class Pen:
    def __init__(self,D,H):
        self.D = D 
        self.H = H

    @property
    def price(self):
        price = self.D * self.H * 1000
        return price

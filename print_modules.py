import numpy as np
import modules


class Print_Modules(object):
    def __init__(self, x_in, x_name, p_in):
        self.x = x_in
        self.x_name = x_name
        self.p = p_in
        # merge input dicts
        self.wec, self.wave, self.pen, self.fish = modules.input_merge(self.x, self.x_name, self.p)
        self.J = self.obj
        
    def obj(self):
        return modules.obj(self.x, self.x_name, self.p)
    
    def P_rated(self):
        modules.power(self.wec, self.wave)
        return self.wec.P_rated
    
    def price_breakdown(self):
        modules.power(self.wec, self.wave)
        return self.wec.price, self.pen.price, modules.fish_feed_price(self.pen, self.fish)
    
    
from math import pi
from typing import Dict

class WEC:
    def __init__(self, capture_width: float, 
                capture_width_ratio_dict: Dict[str,float], 
                wave_damping_dict: Dict[str,float], 
                wec_type: str,
                unit_cost: float) -> None:

        self.capture_width = capture_width
        self.capture_width_ratio_dict = capture_width_ratio_dict
        self.wave_damping_dict = wave_damping_dict
        self.wec_type = wec_type
        self.unit_cost = unit_cost

    @property
    def price(self) -> float:
        price = self.capture_width * self.unit_cost
        return price

    @property
    def wave_damping(self) -> float:
        damping = self.capture_width_ratio_dict[self.wec_type]
        return damping

    @property
    def capture_width_ratio(self) -> float:
        capture_width_ratio = self.capture_width_ratio_dict[self.wec_type]
        return capture_width_ratio

class Wave:
    def __init__(self, Hs: float, T: float) -> None:
        self.Hs = Hs
        self.T = T
        self.rho = 1000
        self.g = 9.81
    
    @property
    def power(self) -> float:
        P_wave = 1/32 * 1/pi * self.rho * self.g**2 * self.Hs**2 * self.T
        return P_wave

class Pen:
    def __init__(self, D: float, H: float, unit_cost: float) -> None:
        self.D = D 
        self.H = H
        self.unit_cost = unit_cost

    @property
    def price(self) -> float:
        price = self.D * self.H * self.unit_cost
        return price

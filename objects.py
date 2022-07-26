from math import cos, exp, pi
from typing import Dict
import numpy as np
from scipy.integrate import trapz

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
        
        self.P_gen = []

    @property
    def price(self) -> float:
        #price = self.capture_width * self.unit_cost
        price = self.P_gen * self.unit_cost
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
    def __init__(self, D: float, H: float, SD: float, n: float, spacing: float, 
                 unit_cost: float, loss_rate: float, harvest_weight: float, temp: float, 
                 O2_in: float, O2_min: float, P_f: float, P_p: float, U_min: float, tau: float, 
                 permeability: float, F_f: float, F_p: float, F_c: float, A_f: float, A_p: float,
                 A_c: float, O_f: float, O_p: float, O_c: float, C_f: float, C_p: float, C_c: float) -> None:
        self.D = D 
        self.H = H
        self.SD = SD 
        self.n = n 
        self.spacing = spacing

        self.unit_cost = unit_cost
        self.loss_rate = loss_rate
        self.harvest_weight = harvest_weight
        self.temp = temp
        
        self.O2_in = O2_in
        self.O2_min = O2_min
        self.P_f = P_f
        self.P_p = P_p
        self.U_min = U_min
        self.tau = tau
        self.permeability = permeability

        self.F_f = F_f
        self.F_p = F_p
        self.F_c = F_c

        self.A_f = A_f
        self.A_p = A_p
        self.A_c = A_c
        
        self.O_f = O_f
        self.O_p = O_p
        self.O_c = O_c
        
        self.C_f = C_f
        self.C_p = C_p
        self.C_c = C_c
        
        self.fish_yield = []
    
    @property
    def price(self) -> float:
        price = self.fish_yield * self.unit_cost
        #price = self.D * self.H * self.unit_cost
        return price

    @property 
    def volume(self) -> float:
        volume = pi * self.D**2 / 4 * self.H
        return volume

    @property 
    def DO2(self) -> float:
        # specific energy content of feed
        delta = self.F_f * self.C_f + self.F_p * self.C_p + self.F_c * self.C_c

        # specific energy content of fish
        C_f_star = 0.85 * self.C_p * self.P_p + self.C_f * self.P_f

        # fraction of food energy contributed by progeins, fat, and carbs
        E_p = self.F_p * self.C_p / delta
        E_f = self.F_f * self.C_f / delta
        E_c = self.F_c * self.C_c / delta

        # metabolizable energy content of food
        FL = (1-self.A_p) * E_p + (1-self.A_f) * E_f + (1-self.A_c) * E_c
        BC = 0.3 * self.A_p * E_p + 0.05 * (self.A_f * E_f + self.A_c * E_c)
        eps = 1 - FL - BC
        eps_star = eps - 0.15 * self.F_p * self.C_p * self.A_p / delta

        # Water temperature as a function of time
        time = np.linspace(0,51,1) # time vector [weeks]
        T = 52                  # period [weeks]
        w = 2*pi/T              # frequency [1/weeks]
        phi = 2*pi/3            # phase offset [-]
        T_max = 23
        T_min = 4
        T_bar = (T_max+T_min)/2
        T_amp = (T_max-T_min)/2
        Temp = self.temp #T_bar + T_amp * cos(w * time + phi)

        # Fish growth as a function of time
        a = 0.038
        W_0 = 0
        #print(Temp)
        #integral = trapz( exp(Temp*self.tau), x=time )
        integral = Temp # fixme
        W = (W_0**(1/3) + a/3 * integral)**3

        # Growth rate as a function of time
        b = 2/3
        W_dot = a * W**b * exp(Temp * self.tau)

        # Rate of energy ingested by fish, cal/day
        alpha = 11
        gamma = 0.8
        Q_r = 1/eps_star * (alpha * W**gamma + a * C_f_star * W**b) * exp(time*self.tau)

        # Respiratory oxygen demand with respect to protein, fat, and carb consumption of fish
        DO2_p = (self.F_p * self.A_p * Q_r / delta - self.P_p * W_dot) * self.O_p
        DO2_f = (self.F_f * self.A_f * Q_r / delta - self.P_f * W_dot) * self.O_f
        DO2_c = self.F_c * self.A_c * Q_r / delta * self.O_c

        # Total respiratory oxygen demand of fish per day
        DO2 = DO2_p + DO2_f + DO2_c

        return DO2

    @property
    def carrying_capacity(self) -> float:
        length = self.n * self.D + self.spacing * (self.n-1)
        carrying_capacity = (self.O2_in - self.O2_min) * length * self.H * self.permeability * self.U_min / self.DO2
        #print('carrying capacity: ', carrying_capacity)
        return min([carrying_capacity])


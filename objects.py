from math import cos, exp, pi
from typing import Dict
import numpy as np
from scipy.integrate import trapz
from scipy.integrate import quad
import math
from matplotlib import pyplot as plt

    
class WEC:
    def __init__(self, capture_width: Dict[str,float], 
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
    def annual_energy(self) -> float:
        beta_wec = 0.95 * 0.98                      #For RM3 (device availability * transmission efficiency)
        AEP = self.P_gen * 0.001 * 8766 * beta_wec  #Annual Energy Production [kWh]
        return AEP

    @property
    def price(self) -> float:
        price = self.annual_energy * self.unit_cost
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

class Fish:
    def __init__(self, F_f, F_p, F_c, A_f, A_p, A_c, 
                 O_f, O_p, O_c, C_f, C_p, C_c, 
                 P_f, P_p, tau, loss_rate, harvest_weight, 
                 O2_min, U_min, U_max, temp_min, temp_max, 
                 salinity_min, salinity_max) :
               
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
        
        self.P_f = P_f
        self.P_p = P_p
        self.tau = tau
        self.loss_rate = loss_rate
        self.harvest_weight = harvest_weight
        
        self.O2_min = O2_min
        self.U_min = U_min
        self.U_max = U_max
        self.temp_min = temp_min
        self.temp_max = temp_max
        self.salinity_min = salinity_min
        self.salinity_max = salinity_max
        
        self.time_i = []
        self.W_i = []
        self.W_dot_i = []
        self.DO2_p_i = []
        self.DO2_f_i = []
        self.DO2_c_i = []
        self.DO2_i = []
        self.carrying_capacity_i = []
    
    def integrand(self, t):
        return math.exp(self.temp * self.tau)
    
    def DO2(self, temp) -> float:
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
        time_i = np.linspace(1, 365, 365) 
        self.time_i = time_i;
        
        self.temp = temp
        Temp = temp  

        # Fish growth as a function of time
        a = 0.038
        W_0 = 0 #[g]
        W_i = np.zeros(len(time_i))
        W_dot_i = np.zeros(len(time_i))
        DO2_p_i = np.zeros(len(time_i))
        DO2_f_i = np.zeros(len(time_i))
        DO2_c_i = np.zeros(len(time_i))
        DO2_i = np.zeros(len(time_i))
        
        for i in range(len(time_i)):
            time = time_i[i]
            
            integral = quad(self.integrand, 0, time)  
            W = (W_0**(1/3) + a/3 * integral[0])**3
            W_i[i] = W
            
            # Growth rate as a function of time
            b = 2/3
            W_dot = a * W**b * exp(Temp * self.tau)
            W_dot_i[i] = W_dot
            
            # Rate of energy ingested by fish, cal/day
            alpha = 11
            gamma = 0.8
            Q_r = 1/eps_star * (alpha * W**gamma + a * C_f_star * W**b) * exp(Temp*self.tau)
            
            # Respiratory oxygen demand with respect to protein, fat, and carb consumption of fish
            DO2_p = (self.F_p * self.A_p * Q_r / delta - self.P_p * W_dot) * self.O_p
            DO2_f = (self.F_f * self.A_f * Q_r / delta - self.P_f * W_dot) * self.O_f
            DO2_c = self.F_c * self.A_c * Q_r / delta * self.O_c

            # Total respiratory oxygen demand of fish per day [g/day]
            DO2 = DO2_p + DO2_f + DO2_c
            
            DO2_p_i[i] = DO2_p
            DO2_f_i[i] = DO2_f
            DO2_c_i[i] = DO2_c
            DO2_i[i] = DO2 # [g/day] for a fish
        
        self.W_i = W_i
        self.W_dot_i = W_dot_i
        self.DO2_p_i = DO2_p_i
        self.DO2_f_i = DO2_f_i
        self.DO2_c_i = DO2_c_i
        self.DO2_i = DO2_i
        
        #print("DO2_i=", DO2_i)

        W_i_50g = next(x[0] for x in enumerate(self.W_i) if x[1] > 50)
        W_i_1kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 1000) # index of weight fish = 1 kg
        OCR = sum(self.DO2_i[W_i_50g:W_i_1kg]) #np.cumsum(self.DO2_i)[W_i_1kg]
        #print('OCR', OCR)
                
        return OCR

    @property
    def plot_variable(self):
        ax1 = plt.subplot(3,1,1)
        ax1.plot(self.time_i, self.W_i/1000, label='W')
        ax1.set(xlabel='time [day]', ylabel='Fish weight (W [kg])');
        ax1.legend()
        plt.show()
        
        ax2 = plt.subplot(3,1,2)
        ax2.plot(self.time_i, self.DO2_i , label='DO2')
        ax2.set(xlabel='time [day]', ylabel='DO2 [g/day]');
        ax2.legend()
        plt.show()
               
        W_i_50g = next(x[0] for x in enumerate(self.W_i) if x[1] > 50)
        W_i_1kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 1000)
        W_i_2kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 2000)
        W_i_3kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 3000)
        W_i_4kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 4000)
        
        ax3 = plt.subplot(3,1,3)
        ax3.plot(self.W_i/1000,  np.cumsum(self.DO2_i), 'b' , label='Total DO2')

        ref_DO2 = np.full(shape=(len(self.DO2_i),), fill_value=np.NaN)
        ref_DO2[W_i_1kg] = 445
        ref_DO2[W_i_2kg] = 956
        ref_DO2[W_i_3kg] = 1496
        ref_DO2[W_i_4kg] = 2049
        ax3.plot(self.W_i/1000,  ref_DO2, 'r-o', label='Ref Total DO2')
        ax3.set(xlabel='Fish weight [kg]', ylabel='DO2 [g]');
        ax3.legend()
        plt.show()
             
        print('DO2 for 1kg fish',sum(self.DO2_i[W_i_50g:W_i_1kg]))
        print('DO2 for 2kg fish',sum(self.DO2_i[W_i_50g:W_i_2kg]))
        print('DO2 for 3kg fish',sum(self.DO2_i[W_i_50g:W_i_3kg]))
        print('DO2 for 4kg fish',sum(self.DO2_i[W_i_50g:W_i_4kg]))
        
        print('fish weight after 365 days',self.W_i[-1])        
        
class Pen:
    def __init__(self, D: float, H: float, Depth: float, SD: float, n: float, spacing: float, 
                 unit_cost: float, temp: float, O2_in: float, U: float, salinity: float,
                 permeability: float) -> None:
        self.D = D 
        self.H = H
        self.Depth = Depth
        self.SD = SD 
        self.n = n 
        self.spacing = spacing

        self.unit_cost = unit_cost
        self.temp = temp
        self.O2_in = O2_in
        self.U = U
        self.salinity = salinity
        self.permeability = permeability
        
        self.TPF_O2 = 0
        
    @property 
    def volume(self) -> float:
        volume = pi * self.D**2 / 4 * self.H
        return volume

    @property
    def price(self) -> float:
        price = self.n * self.volume * self.unit_cost
        return price

    @property 
    def power(self) -> float:
        power = self.fish_yield * 0.572  # Annual Energy [kWh] (previously 50000 [W])
        return power

    
    def carrying_capacity(self, fish) -> float:
        length = self.n * self.D   # Based on worst case
        #length = self.n * self.D + self.spacing * (self.n-1) # From reference paper for a row farm
        #length = self.D # for each pen

        #print('length',length)
                    
        OT = (self.O2_in - fish.O2_min) * length * self.Depth * self.permeability * fish.U_min # [g_O2 / s]
               
        #print('OT=' , OT)
        
        self.TPF_O2 = (OT * 3600 * 24 * 365) / fish.DO2(self.temp)  # [kg-fish / year]
        #print('fish.DO2=',fish.DO2(self.temp))
        #print('TPF_O2=',self.TPF_O2)
        
        carrying_capacity = (OT * 3600 * 24) / fish.DO2_i[-1]
        #print('DO2_i=',fish.DO2_i[-1])
        
        self.carrying_capacity_value = carrying_capacity
        #print('carrying_capacity=',carrying_capacity)
        
        return carrying_capacity

    
from math import exp, pi
from typing import Dict
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt
from numba import njit
import copy
from scipy.optimize import minimize



@njit
def cumsum_with_limits_nb(values, uplimit):
    n = len(values)
    res = np.empty(n)
    sum_val = 0
    for i in range(n):
        x = values[i] + sum_val
        if (x <= uplimit):
            res[i] = x
            sum_val = x
        else:
            res[i] = uplimit
            sum_val = uplimit
    return res

# ============================================================================ #
#                       Wave Energy Converter                                  #
# ============================================================================ #

class WEC:
    def __init__(self, Hs, Tp, 
                 capture_width: Dict[str,float], 
                capture_width_ratio_dict: Dict[str,float], 
                wave_damping_dict: Dict[str,float], 
                wec_type: str,
                capacity_factor,
                eta, float_diameter,
                #CapEx_ref, OpEx_ref,
                lifetime, discount_rate, P_wec_rated):
        
        self.wave = Wave(Hs, Tp)

        self.capture_width = capture_width
        self.capture_width_ratio_dict = capture_width_ratio_dict
        self.wave_damping_dict = wave_damping_dict
        self.wec_type = wec_type
        self.capacity_factor = capacity_factor
        self.eta = eta
        self.beta_wec = 0.95 * 0.98  #For RM3 (device availability * transmission efficiency)
        self.float_diameter = float_diameter
        # self.CapEx_ref = CapEx_ref
        # self.OpEx_ref = OpEx_ref
        self.lifetime = lifetime
        self.discount_rate = discount_rate
        self.P_wec_rated = P_wec_rated

        self.P_gen = []
        
    @property
    def AEP(self): # annual energy production for an array of WECs
        AEP = np.sum(self.P_gen)  #Annual Energy Production [kWh]
        return AEP

    @property
    def P_ave(self):
        return self.AEP / 8766 / self.wec_number  #[kW]

    @property
    def CapEx_ref(self):
        array_scale = [1, 10, 50, 100]
        CapEx_array_scale = np.array([61140, 21400, 14490, 13630]) * 300 * 1.29 #'[$/wec_number]'  # inflation rate from 2014 to 2023 = 1.29                                                                           
        f = np.poly1d(np.polyfit(array_scale, CapEx_array_scale, 2))
        return f(self.wec_number)

    @property
    def OpEx_ref(self):
        array_scale = [1, 10, 50, 100]
        OpEx_array_scale = np.array([4070, 1150, 460, 330]) * 300 * 1.29 #'[$/wec_number]'  # inflation rate from 2014 to 2023 = 1.29                                                                           
        f = np.poly1d(np.polyfit(array_scale, OpEx_array_scale, 2))
        return f(self.wec_number)
    
    @property
    def price(self):
        price = self.AEP * self.LCOE
        return price
    
    @property
    def CapEx(self):
        CapEx = self.wec_number * self.CapEx_ref  # capital expense for array
        return CapEx
    
    @property
    def OpEx(self):
        OpEx = self.wec_number * self.OpEx_ref  # operational  expense for array
        return OpEx

    @property
    def NPV_OpEx(self):
        NPV_OpEx = 0
        for i in range(self.lifetime):
            NPV_OpEx += (self.OpEx) / ((1+self.discount_rate)**(i+1))
        return NPV_OpEx
        
    @property
    def cost_NPV(self): #net present value
        cost_NPV = self.CapEx + self.NPV_OpEx
        return cost_NPV

    @property
    def LCOE_base_RM3(self):
        return 0.75 * 1.19   #'[$/kWh]'  # 0.75 [$/kWh] for 100 wecs     inflation rate from 2014 to 2022 = 1.19   
    
    @property
    def LCOE(self):
        return self.LCOE_base_RM3 * (700*1000) / self.AEP_per_unit #  AEP_per_unit for base RM3 = 700 MWh
    
    @property
    def AEP_per_unit(self): # average annual energy production per unit of RM3
        return self.AEP / self.wec_number

    @property
    def wave_damping(self):
        damping = self.capture_width_ratio_dict[self.wec_type]
        return damping

    @property
    def capture_width_ratio(self):
        capture_width_ratio = self.capture_width_ratio_dict[self.wec_type]
        return capture_width_ratio
    
    @property
    def P_mechanical(self):
        P_mechanical = self.wave.P_wave * self.capture_width * self.beta_wec
        return P_mechanical 
    
    @property
    def P_electrical(self):
        P_electrical = np.clip(self.eta * self.P_mechanical, 0, self.wec_number * self.P_wec_rated)
        return P_electrical

    @property
    def wec_number(self):
        number = self.capture_width / (self.float_diameter * self.capture_width_ratio)
        wec_number = np.ceil(number)
        return number

    '''
    @property
    def P_rated(self):
        P_rated = self.P_gen / self.capacity_factor
        return P_rated
    '''

# ============================================================================ #
#                       Wave                                                   #
# ============================================================================ #

class Wave:
    def __init__(self, Hs, Tp) -> None:
        self.Hs = Hs
        self.Tp = Tp  # energy period
        self.rho = 1030
        self.g = 9.81
    
    @property
    def P_wave(self) -> float:
        #P_wave = self.rho * self.g**2 * self.Hs**2 * self.Tp  / (32 * pi)  #for regular wave
        P_wave = self.rho * self.g**2 * self.Hs**2 * self.Tp  / (64 * pi)  * 0.001  # [kw] for irregular wave
        return P_wave

# ============================================================================ #
#                       Fish (Growth and Dissolved Oxygen)                     #
# ============================================================================ #

class Fish:
    def __init__(self, F_f, F_p, F_c, A_f, A_p, A_c, 
                 O_f, O_p, O_c, C_f, C_p, C_c, 
                 P_f, P_p, tau, loss_rate, harvest_weight, 
                 O2_min, U_min, U_max, temp_min, temp_max, 
                 salinity_min, salinity_max, fish_life_cycle,
                 fingerling_weight, fingerling_unit_cost) :
               
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

        self.fish_life_cycle = fish_life_cycle
        self.fingerling_weight = fingerling_weight
        self.fingerling_unit_cost = fingerling_unit_cost
        
        self.time_i = []
        self.W_i = []
        self.W_dot_i = []
        self.DO2_p_i = []
        self.DO2_f_i = []
        self.DO2_c_i = []
        self.DO2_i = []
        self.carrying_capacity_i = []
    
    def integrand(self, t):
        return exp(self.temp * self.tau)
    
    def sigma(self, lower, upper):
        total = 0
        for i in range(lower, upper):
            total += exp(self.temp[i] * self.tau)
        return [total , 0]
    
    def DO2(self, temp):
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
        time_i = np.arange(1, 365, 1)
        self.time_i = time_i
        
        self.temp = temp #[0]

        # Fish growth as a function of time
        a = 0.038
        W_0 = self.fingerling_weight * 1000 #50 #[g]
        W_i = np.zeros(len(time_i))
        W_dot_i = np.zeros(len(time_i))
        DO2_p_i = np.zeros(len(time_i))
        DO2_f_i = np.zeros(len(time_i))
        DO2_c_i = np.zeros(len(time_i))
        DO2_i = np.zeros(len(time_i))
        
        for i in range(len(time_i)):
            time = time_i[i]
            
            integral = self.sigma(0, time) 
            #integral = quad(self.integrand, 0, time)   
            W = (W_0**(1/3) + a/3 * integral[0])**3
            W_i[i] = W
            
            # Growth rate as a function of time
            b = 2/3
            W_dot = a * W**b * exp(self.temp[i] * self.tau)
            #W_dot = a * W**b * exp(self.temp * self.tau)
            W_dot_i[i] = W_dot
            
            # Rate of energy ingested by fish, cal/day
            alpha = 11
            gamma = 0.8
            Q_r = 1/eps_star * (alpha * W**gamma + a * C_f_star * W**b) * exp(self.temp[i]*self.tau)
            #Q_r = 1/eps_star * (alpha * W**gamma + a * C_f_star * W**b) * exp(self.temp*self.tau)
            
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
        
        W_i_50g = next(x[0] for x in enumerate(self.W_i) if x[1] > 50)   # index of weight fish = 50 g
        try:
            W_i_1kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 1000) # index of weight fish = 1 kg
            OCR = sum(self.DO2_i[W_i_50g:W_i_1kg])
        except:
            OCR = sum(self.DO2_i[W_i_50g:]) / self.W_i[-1] * 1000

        return OCR

    def plot_variable(self):
        fig, axes = plt.subplots(4,1, figsize=(12, 11))

        ax0 = axes[0]
        ax0.plot(self.temp)
        ax0.set(xlabel='time [day]', ylabel='Water Temperature [C]');
        #ax0.legend()
        ax0.grid(True)
        ax0.set_xlim(0, None)

        ax1 = axes[1]
        ax1.plot(self.time_i, self.W_i/1000)
        ax1.set(xlabel='time [day]', ylabel='Fish weight (W [kg])');
        #ax1.legend()
        ax1.grid(True)
        ax1.set_xlim(0, None)
        
        ax2 = axes[2]
        ax2.plot(self.time_i, self.DO2_i)
        # ax2.plot(self.time_i, self.DO2_p_i, label='pi')
        # ax2.plot(self.time_i, self.DO2_f_i, label='fi')
        # ax2.plot(self.time_i, self.DO2_c_i, label='ci')
        ax2.set(xlabel='time [day]', ylabel='DO2 [g/day]');
        #ax2.legend()
        ax2.grid(True)
        ax2.set_xlim(0, None)

        ref_DO2 = np.full(shape=(len(self.DO2_i),), fill_value=np.NaN)
        W_i_50g = next(x[0] for x in enumerate(self.W_i) if x[1] > 50)
        W_i_1kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 1000)
        ref_DO2[W_i_1kg] = 445
        print('DO2 for 1kg fish',sum(self.DO2_i[W_i_50g:W_i_1kg]))
        try:
            W_i_2kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 2000)
            ref_DO2[W_i_2kg] = 956
            print('DO2 for 2kg fish',sum(self.DO2_i[W_i_50g:W_i_2kg]))
            try:
                W_i_3kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 3000)
                ref_DO2[W_i_3kg] = 1496
                print('DO2 for 3kg fish',sum(self.DO2_i[W_i_50g:W_i_3kg]))
                try:
                    W_i_4kg = next(x[0] for x in enumerate(self.W_i) if x[1] > 4000)
                    ref_DO2[W_i_4kg] = 2049
                    print('DO2 for 4kg fish',sum(self.DO2_i[W_i_50g:W_i_4kg]))
                except:
                    W_i_4kg = np.nan
            except:
                W_i_3kg = np.nan

        except:
            W_i_2kg = np.nan
        
        print('fish weight after 365 days',self.W_i[-1])        
        
        ax3 = axes[3]
        ax3.plot(self.W_i/1000,  np.cumsum(self.DO2_i), 'b' , label='Total DO2')

        ax3.plot(self.W_i/1000,  ref_DO2, 'r-o', label='Ref Total DO2')
        ax3.set(xlabel='Fish weight [kg]', ylabel='DO2 [g]')
        ax3.legend()
        ax3.grid(True)
        ax3.set_xlim(0, None)

        #plt.subplots_adjust(hspace=0.3)
        plt.tight_layout()
        plt.show()
        return

# ============================================================================ #
#                                Net Pen                                       #
# ============================================================================ #
         
class NetPen:
    def __init__(self, D, H, total_cost=None):
        self.D = D 
        self.H = H
        self.total_cost = total_cost

        self.walk_side = 1 # [m]

    @property 
    def h_cylinder(self):
        h_cylinder = self.H / 1.3 # we assume cone height is 30% of cylinder height
        return h_cylinder
    
    @property 
    def h_cone(self):
        h_cone = self.H - self.h_cylinder
        return h_cone
    
    @property 
    def volume(self):
        #volume = top cylinder_volume + cone_volume
        cylinder_volume = pi * self.D**2 / 4 * self.h_cylinder
        cone_volume = 1/3 * pi * self.D**2 / 4 * self.h_cone 
        volume = cylinder_volume + cone_volume
        return volume
    
    @property 
    def surface(self):
        #surface = top_surface + cone_surface + side_surface
        top_surface = pi * self.D**2 / 4
        cone_surface = pi * self.D/2 * np.sqrt(self.h_cone**2 + self.D**2 / 4)
        side_surface = pi * self.D * self.h_cylinder
        surface = top_surface + cone_surface + side_surface
        return surface
    
    @property 
    def struct_length(self):
        #length = float collar inner_perimeter + float collar outer_perimeter + top net support float_perimeter + sinker tube_perimeter
        float_collar_inner_perimeter = pi * self.D / 2
        float_collar_outer_perimeter = pi * (self.D + 2 * self.walk_side) / 2 
        top_net_support_float_perimeter = pi * self.D * self.walk_side # we assume net_support every 1 meter
        sinker_tube_perimeter = pi * self.D / 2
        struct_length = float_collar_inner_perimeter + float_collar_outer_perimeter + top_net_support_float_perimeter + sinker_tube_perimeter
        return struct_length
    
    @property 
    def netting_cost(self):
        netting_cost = 0.35 * self.total_cost
        return netting_cost

    @property 
    def struct_cost(self):
        struct_cost = 0.65 * self.total_cost
        return struct_cost

        
# ============================================================================ #
#         Aquacalture Farm (Net Pen, Feed Barge, Carrying Capacity)            #
# ============================================================================ #
         
class Aqua:
    def __init__(self, 
                 F_f, F_p, F_c, A_f, A_p, A_c, 
                 O_f, O_p, O_c, C_f, C_p, C_c, 
                 P_f, P_p, tau, loss_rate, harvest_weight, 
                 O2_min, U_min, U_max, temp_min, temp_max, 
                 salinity_min, salinity_max, fish_life_cycle,
                 fingerling_weight, fingerling_unit_cost,

                 D, #H, 
                 Depth, SD, pen_number, spacing, 
                 temp, O2_in, U, salinity, permeability, bathymetry,
                 pos_lat, pos_long, 
                 pen_netting_CapEx_ref, pen_struct_CapEx_ref, feedbarge_CapEx_ref, feedbarge_OpEx_ref,
                 lifetime, discount_rate,
                 FCR, feed_unit_cost, feedbarge_unit_capacity, feedbarge_unit_feedlines,
                 summer_feedbarge_power, summer_lighting_power_per_kg, summer_equipment_power_per_kg,
                 winter_feedbarge_power, winter_lighting_power_per_kg, winter_equipment_power_per_kg,
                 
                 fuel_consump_rate, fuel_cost, 
                 captain_salary, crew_salary, crew_num, 
                 t_feed, velocity, number_travel, distance
                 ):
        
        self.netpen_ref_model()

        self.fish = Fish(F_f, F_p, F_c, A_f, A_p, A_c, 
                 O_f, O_p, O_c, C_f, C_p, C_c, 
                 P_f, P_p, tau, loss_rate, harvest_weight, 
                 O2_min, U_min, U_max, temp_min, temp_max, 
                 salinity_min, salinity_max, fish_life_cycle,
                 fingerling_weight, fingerling_unit_cost)
        
        H = self.netpen_geometry(D)
        self.netpen = NetPen(D, H)
        
        self.vessel = Vessel(fuel_consump_rate, fuel_cost, 
                 captain_salary, crew_salary, crew_num, 
                 t_feed, velocity, number_travel, distance,
                 lifetime, discount_rate)
        
       
        self.Depth = Depth
        self.SD = SD 
        self.pen_number = pen_number
        self.spacing = spacing
        
        self.waterdepth_underpen_min = 20 # [m]
        self.waterdepth_underpen_max = 75 # [m]

        self.temp = temp
        self.O2_in = O2_in
        self.U = U
        self.salinity = salinity
        self.permeability = permeability
        
        self.bathymetry = bathymetry
        self.waterdepth_underpen_min = 20 # [m]
        self.waterdepth_underpen_max = 75 # [m]

        self.pos_lat = pos_lat
        self.pos_long = pos_long

        self.lifetime = lifetime
        self.discount_rate = discount_rate

        self.feedbarge_CapEx_ref = feedbarge_CapEx_ref
        self.feedbarge_OpEx_ref = feedbarge_OpEx_ref
        self.feedbarge_unit_capacity = feedbarge_unit_capacity
        self.feedbarge_unit_feedlines = feedbarge_unit_feedlines

        self.FCR = FCR
        self.feed_unit_cost = feed_unit_cost
        
        self.TPF_O2 = 0

        self.summer_feedbarge_power = summer_feedbarge_power
        self.summer_lighting_power_per_kg = summer_lighting_power_per_kg
        self.summer_equipment_power_per_kg = summer_equipment_power_per_kg
        self.winter_feedbarge_power = winter_feedbarge_power
        self.winter_lighting_power_per_kg = winter_lighting_power_per_kg
        self.winter_equipment_power_per_kg = winter_equipment_power_per_kg

    def netpen_ref_model(self):
        ## data from INNOVASEA
        pen5200  = NetPen(D=24, H=13, total_cost=365000)
        pen7600  = NetPen(D=27, H=16, total_cost=405000)
        pen14000 = NetPen(D=33, H=20, total_cost=595000)

        x1 = np.array([pen5200.surface,      pen7600.surface,      pen14000.surface])
        y1 = np.array([pen5200.netting_cost, pen7600.netting_cost, pen14000.netting_cost])
        netting_cost_coeff = np.polyfit(x1, y1, 1)
        self.netting_cost = np.poly1d(netting_cost_coeff)

        x2 = np.array([pen5200.D,           pen7600.D,           pen14000.D])
        y2 = np.array([pen5200.struct_cost, pen7600.struct_cost, pen14000.struct_cost])
        struct_cost_coeff = np.polyfit(x2, y2, 1)
        self.struct_cost = np.poly1d(struct_cost_coeff)

        x3 = np.array([pen5200.D,      pen7600.D,      pen14000.D])
        y3 = np.array([pen5200.H,      pen7600.H,      pen14000.H])
        netpen_geometry_coeff = np.polyfit(x3, y3, 1)
        self.netpen_geometry = np.poly1d(netpen_geometry_coeff)
    
    @property
    def fish_feed_harvest_week(self):
        fish_feed_harvest_week = self.fish_yield * .02 * 7  # require a daily fish feed of 2% of fish yield for a week
        return fish_feed_harvest_week

    @property 
    def feedbarge_number(self):
        if ((self.fish_feed_harvest_week / self.pen_number) < (self.feedbarge_unit_capacity / self.feedbarge_unit_feedlines)):
            feedbarge_number = self.pen_number / self.feedbarge_unit_feedlines
        else:
            #feedbarge_number = np.ceil(self.fish_feed_harvest_week / self.feedbarge_unit_capacity)
            feedbarge_number = self.fish_feed_harvest_week / self.feedbarge_unit_capacity
        return feedbarge_number

    @property
    def CapEx_feedbarge(self):
        CapEx_feedbarge = self.feedbarge_number * self.feedbarge_CapEx_ref
        return CapEx_feedbarge

    @property
    def CapEx_netpen_struct(self):
        CapEx_netpen_struct = self.struct_cost(self.netpen.D)
        return CapEx_netpen_struct

    @property
    def CapEx_netpen_netting(self):
        CapEx_netpen_netting = self.netting_cost(self.netpen.surface)
        return CapEx_netpen_netting
    
    @property
    def Capex_netpen_mooring(self):
        Capex_netpen_mooring = 0.25 * (self.CapEx_netpen_struct + self.CapEx_netpen_netting)
        return Capex_netpen_mooring

    @property
    def CapEx_netpen(self):
        CapEx_per_pen = self.CapEx_netpen_struct + self.CapEx_netpen_netting + self.Capex_netpen_mooring
        #CapEx_per_pen = (self.netpen.volume * 100)
        CapEx_netpen = self.pen_number * CapEx_per_pen
        return CapEx_netpen

    @property
    def CapEx(self):
        CapEx = self.CapEx_netpen + self.CapEx_feedbarge
        return CapEx

    @property
    def NPV_OpEx(self):
        NPV_OpEx = self.NPV_OpEx_fish_feed_price + self.NPV_OpEx_fingerling_price_annual + self.vessel.NPV_OpEx
        return NPV_OpEx
    
    @property
    def NPV_OpEx_fish_feed_price(self):
        NPV_OpEx_fish_feed_price = 0
        for i in range(self.lifetime):
            NPV_OpEx_fish_feed_price += (self.fish_feed_price_annual) / ((1+self.discount_rate)**(i+1))
        return NPV_OpEx_fish_feed_price

    @property
    def NPV_OpEx_fingerling_price_annual(self):
        NPV_OpEx_fingerling_price_annual = 0
        for i in range(self.lifetime):
            NPV_OpEx_fingerling_price_annual += (self.fingerling_price_annual) / ((1+self.discount_rate)**(i+1))
        return NPV_OpEx_fingerling_price_annual
        
    @property
    def cost_NPV(self): #net present value
        cost_NPV = self.CapEx + self.NPV_OpEx
        return cost_NPV
    
    @property 
    def power_summer(self):
        power_summer = (self.summer_feedbarge_power * self.feedbarge_number) + (self.summer_lighting_power_per_kg * self.fish_yield) + (self.summer_equipment_power_per_kg * self.fish_yield)
        return power_summer

    @property 
    def power_winter(self):
        power_winter = (self.winter_feedbarge_power * self.feedbarge_number) + (self.winter_lighting_power_per_kg * self.fish_yield) + (self.winter_equipment_power_per_kg * self.fish_yield)
        return power_winter
    
    @property 
    def power(self):
        summer_weight = np.zeros(365)
        winter_weight = np.zeros(365)
        power = []
        season_length = 91
        for i in range(365):
            if (i < season_length): # winter
                summer_weight[i] = 0
            elif (i < 2*season_length): # spring
                summer_weight[i] = min(summer_weight[i-1] + (1 / season_length), 1)
            elif (i < 3*season_length): # summer
                summer_weight[i] = 1
            else: #fall
                summer_weight[i] = max(summer_weight[i-1] - (1 / season_length) , 0)

            winter_weight[i] = 1 - summer_weight[i]
            
            daily_power = summer_weight[i] * self.power_summer + winter_weight[i] * self.power_winter

            power = np.append(power, daily_power)
        #noise= np.random.rand(8760)
        #power = np.ones(8760) * (self.fish_yield / 1e5)
        #power[2800:3900] = power[2800:3900] * 1.1
        #power[6900:8000] = power[6900:8000] * 0.9
        return power

    @property
    def fingerling_price(self):
        fingerling_price = self.fish_yield / self.fish.harvest_weight / (1-self.fish.loss_rate) * self.fish.fingerling_unit_cost 
        return fingerling_price
    
    @property
    def fingerling_price_annual(self):
        fingerling_price_annual = self.fingerling_price * 365/self.fish.fish_life_cycle
        return fingerling_price_annual
    
    @property
    def biomass(self):
        biomass = self.pen_number * self.SD * self.netpen.volume  # [kg]
        return biomass
    
    @property
    def fish_yield(self):
        fish_yield = self.pen_number * self.SD * self.netpen.volume  # [kg]
        return fish_yield

    @property
    def fish_feed_price(self):
        fish_feed_price = self.fish_yield * self.FCR * self.feed_unit_cost
        return fish_feed_price
    
    @property
    def fish_feed_price_annual(self):
        fish_feed_price_annual = self.fish_feed_price * 365/self.fish.fish_life_cycle
        return fish_feed_price_annual

    @property
    def carrying_capacity(self):
        length = self.pen_number * self.netpen.D
                    
        OT = (self.O2_in - self.fish.O2_min) * length * self.netpen.H * self.permeability * self.fish.U_min # [g_O2 / s]
        self.TPF_O2 = (OT * 3600 * 24 * 365) / self.fish.DO2(self.temp)  # [kg-fish / year]

        carrying_capacity = (OT * 3600 * 24) / self.fish.DO2_i[-1] 
        
        return carrying_capacity

# ============================================================================ #
#                              Vessel Travel                                   #
# ============================================================================ #

class Vessel:
    def __init__(self, fuel_consump_rate, fuel_cost, 
                 captain_salary, crew_salary, crew_num, 
                 t_feed, velocity, number_travel, distance,
                 lifetime, discount_rate):
        self.fuel_consump_rate = fuel_consump_rate
        self.fuel_cost = fuel_cost
        self.captain_salary = captain_salary
        self.crew_salary = crew_salary
        self.crew_num = crew_num
        self.t_feed = t_feed
        self.velocity = velocity
        self.number_travel = number_travel
        self.distance = distance
        self.t_travel = 2 * self.distance / self.velocity  # back and forth travel between port and deployment location
        self.lifetime = lifetime
        self.discount_rate = discount_rate

    @property 
    def OpEx(self):
        fuel_price = self.t_travel * self.fuel_cost * self.fuel_cost
        laber_salary = (self.t_travel) * (self.captain_salary + self.crew_num * self.crew_salary)
        OpEx = self.number_travel  * (fuel_price + laber_salary) # annual price
        return OpEx
    
    @property
    def NPV_OpEx(self): #net present value
        NPV_OpEx = 0
        for i in range(self.lifetime):
            NPV_OpEx += (self.OpEx) / ((1+self.discount_rate)**(i+1))
        return NPV_OpEx

# ============================================================================ #
#                       Diesel Generator                                       #
# ============================================================================ #

class DieselGen:
    def __init__(self, fuel_consump_rate, fuel_cost, eta, load_level, CapEx_ref, OpEx_ref, lifetime, discount_rate):
        self.fuel_consump_rate = fuel_consump_rate
        self.fuel_cost = fuel_cost
        self.eta = eta
        self.load_level = load_level
        self.CapEx_ref = CapEx_ref
        self.OpEx_ref = OpEx_ref
        self.lifetime = lifetime
        self.discount_rate = discount_rate
        self.P_rated = []
    
    def power(self, power_out):
        self.P_rated = power_out / self.eta / self.load_level
        return
    
    @property
    def CapEx(self):
        CapEx = self.P_rated  * 8760 * self.CapEx_ref # capital expense 
        return CapEx
    
    @property
    def OpEx(self):
        OpEx = (self.fuel_cost * self.fuel_consump_rate) * 8760 + self.P_rated  * 8760 * self.OpEx_ref  # operational expense
        return OpEx

    @property
    def NPV_OpEx(self): #net present value
        NPV_OpEx = 0
        for i in range(self.lifetime):
            NPV_OpEx += (self.OpEx) / ((1+self.discount_rate)**(i+1))
        return NPV_OpEx
        
    @property
    def cost_NPV(self): #net present value
        cost_NPV = self.CapEx + self.NPV_OpEx
        return cost_NPV

# ============================================================================ #
#                              Energy Storage                                  #
# ============================================================================ #

class ES:
    def __init__(self, eta, CapEx_ref, OpEx_ref, lifetime, discount_rate, soc_uplimit, soc_downlimit):
        self.eta = eta
        self.CapEx_ref = CapEx_ref
        self.OpEx_ref = OpEx_ref
        self.lifetime = lifetime
        self.discount_rate = discount_rate
        self.soc_uplimit = soc_uplimit
        self.soc_downlimit = soc_downlimit
        self.P_diff = []
        self.size = []
        self.total_size = []
        self.power = []

    #Hourly Power for Energy Storage
    @property 
    def P_stored_hour(self):
        P_stored_hour = np.zeros(len(self.P_diff))
        for i in range(len(self.P_diff)):
            if self.P_diff[i]>0:
                P_stored_hour[i] = self.P_diff[i] * self.eta
            else:
                P_stored_hour[i] = self.P_diff[i] / self.eta
        P_stored_hour[0] += self.size # assuming the battery is fully charged at the start time
        return P_stored_hour
    
    #Cummulative Power at Energy Storage
    @property 
    def P_stored_cum(self):
        P_stored_cum = cumsum_with_limits_nb(self.P_stored_hour, self.size * self.soc_uplimit)
        return P_stored_cum
    
    def sizing_func(self, P_diff):
        self.P_diff = P_diff
        self.size  = 0
        self.total_size = abs(np.min(self.P_stored_cum)) / (self.soc_uplimit - self.soc_downlimit)
    
        self.size = self.total_size
        self.power = copy.deepcopy(self.P_stored_cum)
    
    @property
    def soc(self):
        soc = self.power / self.total_size
        return soc

    @property
    def CapEx(self):
        CapEx = self.total_size * self.CapEx_ref  # capital expense
        return CapEx
    
    @property
    def OpEx(self):
        OpEx = self.total_size * self.OpEx_ref  # operational expense 
        return OpEx
    
    @property
    def cost_NPV(self): #net present value
        cost_NPV = self.CapEx
        for i in range(self.lifetime):
            cost_NPV += (self.OpEx) / ((1+self.discount_rate)**(i+1))
        return cost_NPV #/ 1000000 #[M$]
    
    def plot_es(self):
        plt.plot(self.power)
        plt.axhline(y=self.total_size, color='r', linestyle='-')
        plt.axhline(y=self.total_size * self.soc_uplimit, color='green', linestyle='dotted')
        plt.axhline(y=self.total_size * self.soc_downlimit, color='green', linestyle='dotted')
        plt.axhline(y=0, color='r', linestyle='-')
        plt.ylabel('P_stored [W]')
        plt.xlabel('Time [h]')
        plt.grid()
        plt.show()
    
    def plot_es_soc(self):
        plt.plot(self.soc * 100)
        plt.axhline(y=100, color='r', linestyle='-')
        plt.axhline(y=0, color='r', linestyle='-')
        plt.ylabel('SOC [%]')
        plt.xlabel('Time [h]')
        plt.grid()
        plt.show()
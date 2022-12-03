import numpy as np
import pandas as pd
from modules import Aqua_Obj

def print_bold(str):
    print('\033[1m' + str + '\033[0;0m')
    return

def print_objective(title, aqua_obj):
    print_bold(title+" objective function terms:")
    print(' '*2, "cost_per_yield", "{:10.3f}".format(aqua_obj.cost_per_yield), '[$/kg]')
    print(' '*2, "price         ", "{:10.3f}".format(aqua_obj.price), '[$]')
    print(' '*2, "fish_yield    ", "{:10.3f}".format(aqua_obj.fish_yield), '[kg]')
    print("-"*40)

def print_P_rated(title, aqua_obj):
    print_bold(title+" WEC rated power:")
    print(' '*2, "P_rated     ", "{:10.3f}".format(aqua_obj.wec.P_rated), '[kW]')
    print("-"*40)

def print_travel_distance(title, aqua_obj):
    print_bold(title+" Port to deployment location distance:")
    print(' '*2, "distance     ", "{:10.3f}".format(aqua_obj.env.distance), '[km]')
    print("-"*40)

def print_price_breakdown(title, aqua_obj):
    print_bold(title+" price break down:")
    print(' '*2, "wec price           ", "{:10.3f}".format(aqua_obj.wec.price), '[$]')
    print(' '*2, "pen price           ", "{:10.3f}".format(aqua_obj.pen.price), '[$]')
    print(' '*2, "fish feed price     ", "{:10.3f}".format(aqua_obj.fish_feed_price), '[$]')
    print(' '*2, "energy st price     ", "{:10.3f}".format(aqua_obj.es.price), '[$]')
    print(' '*2, "vessel travel price ", "{:10.3f}".format(aqua_obj.vessel.price), '[$]')
    print("-"*40)

def print_carrying_capacity(title, aqua_obj):
    print_bold(title+" carrying capacity:")
    print(' '*2, "TPF_O2              ", "{:10.3f}".format(aqua_obj.pen.TPF_O2), '[kg fish/year]')
    print(' '*2, "Carrying Capacity   ", "{:10.3f}".format(aqua_obj.carrying_capacity), '[kg fish]')
    print("-"*40)
    
def print_ineq_cons(title,aqua_obj):
    print_bold(title+" constraints:")
    print(' '*2, "P_gen_cons          ", "{:10.3f}".format(aqua_obj.P_gen_cons), '[kWh]')
    print(' '*2, "fish_yield_cons     ", "{:10.3f}".format(aqua_obj.fish_yield_cons), '[kg]')
    print(' '*2, "env_Umin_cons       ", "{:10.3f}".format(aqua_obj.env_Umin_cons), '[m/s]')
    print(' '*2, "env_Umax_cons       ", "{:10.3f}".format(aqua_obj.env_Umax_cons), '[m/s]')
    print(' '*2, "env_tempmin_cons    ", "{:10.3f}".format(aqua_obj.env_tempmin_cons), '[C]')
    print(' '*2, "env_tempmax_cons    ", "{:10.3f}".format(aqua_obj.env_tempmax_cons), '[C]')
    print(' '*2, "env_salinitymin_cons", "{:10.3f}".format(aqua_obj.env_salinitymin_cons), '[PSU]')
    print(' '*2, "env_salinitymax_cons", "{:10.3f}".format(aqua_obj.env_salinitymax_cons), '[PSU]')
    print(' '*2, "env_O2_min_cons     ", "{:10.3f}".format(aqua_obj.env_O2_min_cons), '[mg/l]')
    print("-"*40)
            
def init_result(x0_init, x_name, p_init):
    aqua_init_obj = Aqua_Obj(x0_init, x_name, p_init) 
    print_objective("Initial",aqua_init_obj)
    print_P_rated("Initial",aqua_init_obj)
    print_travel_distance("Initial",aqua_init_obj)
    print_price_breakdown("Initial",aqua_init_obj)
    print_ineq_cons("Initial",aqua_init_obj)
    print_carrying_capacity("Initial",aqua_init_obj)
    print('+'*40)

def optimize_result(x_name, x_list, x_unit, res_opt, p):
    print(res_opt.success)
    aqua_opt_obj = Aqua_Obj(res_opt.x, x_name, p) 
    col_width = len(max(x_list, key=len))
    
    print_bold("optimal design variable:")
    for i in range(len(x_list)):
        print(' '*2, x_list[i], ' '*(col_width - len(x_list[i])) , "{:10.3f}".format(res_opt.x[i]), x_unit[i])
    
    print(' '*2, "num_pens      ", "{:10.3f}".format(aqua_opt_obj.pen.n), '[-]')
    print("-"*40)
    
    print_objective("optimal",aqua_opt_obj)
    print_P_rated("optimal",aqua_opt_obj)
    print_travel_distance("optimal",aqua_opt_obj)
    print_price_breakdown("optimal",aqua_opt_obj)
    print_ineq_cons("optimal",aqua_opt_obj)
    print_carrying_capacity("optimal",aqua_opt_obj)
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
    print(' '*2, "fish_yield    ", "{:10.3f}".format(aqua_obj.pen.fish_yield), '[kg]')
    print("-"*40)

def print_P_rated(title, aqua_obj):
    print_bold(title+" WEC rated power:")
    print(' '*2, "wave power  ", "{:10.3f}".format(aqua_obj.wave_in.P_wave/1000), '[kW]')
    print(' '*2, "P_rated     ", "{:10.3f}".format(aqua_obj.wec.P_rated/1000), '[kW]')
    print("-"*40)

def print_carrying_capacity(title, aqua_obj):
    print_bold(title+" carrying capacity:")
    print(' '*2, "TPF_O2              ", "{:10.3f}".format(aqua_obj.pen.TPF_O2), '[kg fish/year]')
    print(' '*2, "Carrying Capacity   ", "{:10.3f}".format(aqua_obj.carrying_capacity), '[kg fish]')
    print("-"*40)

def print_price_breakdown(title, aqua_obj):
    print_bold(title+" price break down:")
    print(' '*2, "wec price           ", "{:10.3f}".format(aqua_obj.wec.price), '[$]')
    print(' '*2, "pen price           ", "{:10.3f}".format(aqua_obj.pen.price), '[$]')
    print(' '*2, "fish feed price     ", "{:10.3f}".format(aqua_obj.fish_feed_price), '[$]')
    print(' '*2, "vessel travel price ", "{:10.3f}".format(aqua_obj.vessel.price), '[$]')
    print("-"*40)
    
def print_ineq_cons(title,aqua_obj):
    print_bold(title+" constraints:")
    print(' '*2, "normalized P_gen_cons          ", "{:10.3f}".format(aqua_obj.P_gen_cons), '[-]')
    print(' '*2, "normalized fish_yield_cons     ", "{:10.3f}".format(aqua_obj.fish_yield_cons), '[-]')
    print(' '*2, "normalized pen_ratio_low_cons  ", "{:10.3f}".format(aqua_obj.pen_ratio_low_cons), '[-]')
    print(' '*2, "normalized pen_ratio_up_cons   ", "{:10.3f}".format(aqua_obj.pen_ratio_up_cons), '[-]')
    print("-"*40)

def print_point_validation(title,aqua_obj):
    print_bold(title+" point validation:")
    print(' '*2, "conditions          ", aqua_obj.gis_data["ok-conditions"].bool(), '[-]')
    print(' '*2, "scope               ", aqua_obj.gis_data["ok-scope"].bool(), '[-]')
    print(' '*2, "conflicts           ", aqua_obj.gis_data["ok-conflicts"].bool(), '[-]')
    
            
def init_result(x0_init, x_name, p_init):
    aqua_init_obj = Aqua_Obj(x0_init, x_name, p_init) 
    print('+'*40) 
    if (aqua_init_obj.valid_point):
        print_objective("Initial",aqua_init_obj)
        print_P_rated("Initial",aqua_init_obj)
        print_price_breakdown("Initial",aqua_init_obj)
        print_ineq_cons("Initial",aqua_init_obj)
        print_carrying_capacity("Initial",aqua_init_obj)
        print_point_validation("Initial",aqua_init_obj)
    else:
        print('invalid init point')
    print('+'*40) 

def optimize_result(x_name, x_list, x_unit, res_opt, p):
    aqua_opt_obj = Aqua_Obj(res_opt.x, x_name, p) 
    if (aqua_opt_obj.valid_point):
        print('optimization success: ',res_opt.success)
        col_width = len(max(x_list, key=len))
        
        print_bold("optimal design variable:")
        for i in range(len(x_list)):
            print(' '*2, x_list[i], ' '*(col_width - len(x_list[i])) , "{:10.3f}".format(res_opt.x[i]), x_unit[i])
        print("-"*40)
        
        print_objective("optimal",aqua_opt_obj)
        print_P_rated("optimal",aqua_opt_obj)
        print_price_breakdown("optimal",aqua_opt_obj)
        print_ineq_cons("optimal",aqua_opt_obj)
        print_carrying_capacity("optimal",aqua_opt_obj)
        print_point_validation("optimal",aqua_opt_obj)
    else:
        print('invalid optimized point')
    print('+'*40) 

def bruteforce_result(x_name, x_list, x_unit, res_opt, p):
    if (res_opt.valid_point):
        col_width = len(max(x_list, key=len))
        
        print_bold("brute force design variable:")
        for i in range(len(x_list)):
            print(' '*2, x_list[i], ' '*(col_width - len(x_list[i])) , "{:10.3f}".format(res_opt.x0[i]), x_unit[i])
        print("-"*40)
        
        print_objective("brute force",res_opt)
        print_P_rated("brute force",res_opt)
        print_price_breakdown("brute force",res_opt)
        print_ineq_cons("brute force",res_opt)
        print_carrying_capacity("brute force",res_opt)
        print_point_validation("brute force",res_opt)
    else:
        print('invalid optimized point')
    print('+'*40) 

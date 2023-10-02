import numpy as np
import pandas as pd
from modules import Aqua_Obj

def print_bold(str):
    print('\033[1m' + str + '\033[0;0m')
    return

def print_objective(title, aqua_obj):
    print_bold(title+" objective function terms:")
    print(' '*2, "Objective_func", "{:10.3f}".format(aqua_obj.obj_func))
    print(' '*2, "cost per yield", "{:10.3f}".format(aqua_obj.cost_per_yield), '[Million $ / Tonne]')
    print(' '*2, "NPV           ", "{:10.3f}".format(aqua_obj.cost_NPV / 1000000), '[Million $]')
    print(' '*2, "fish yield    ", "{:10.3f}".format(aqua_obj.pen.fish_yield / 1000000), '[Tonne]')
    print("-"*40)

def print_P_rated(title, aqua_obj):
    print_bold(title+" wave energy converter:")
    print(' '*2, "wec number    ", "{:10.3f}".format(aqua_obj.wec.wec_number), '[-]')
    #print(' '*2, "wave power    ", "{:10.3f}".format(aqua_obj.wave_in.P_wave/1000), '[kW/m]')
    print(' '*2, "wec_P_ave     ", "{:10.3f}".format(aqua_obj.wec.P_ave), '[kW]')
    print(' '*2, "wec AEP       ", "{:10.3f}".format(aqua_obj.wec.AEP), '[kWh]')
    print("-"*40)

def print_carrying_capacity(title, aqua_obj):
    print_bold(title+" carrying capacity:")
    print(' '*2, "TPF_O2              ", "{:10.3f}".format(aqua_obj.pen.TPF_O2), '[kg fish/year]')
    print(' '*2, "Carrying Capacity   ", "{:10.3f}".format(aqua_obj.pen.carrying_capacity), '[kg fish]')
    print("-"*40)

def print_price_breakdown(title, aqua_obj):
    print_bold(title+" price break down:")
    print(' '*2, "wec NPV             ", "{:10.3f}".format(aqua_obj.wec.cost_NPV / 1000000), '[Million $]', "{:2.1f}".format(100*aqua_obj.wec.cost_NPV/aqua_obj.cost_NPV), '[%]')
    print(' '*2, "wec CapEx           ", "{:10.3f}".format(aqua_obj.wec.CapEx / 1000000), '[Million $]')
    print(' '*2, "wec OpEx            ", "{:10.3f}".format(aqua_obj.wec.OpEx / 1000000), '[Million $]')

    print(' '*2, "Aqua NPV            ", "{:10.3f}".format(aqua_obj.pen.cost_NPV / 1000000), '[Million $]', "{:2.1f}".format(100*aqua_obj.pen.cost_NPV/aqua_obj.cost_NPV), '[%]')
    print(' '*2, "Net Pen CapEx       ", "{:10.3f}".format(aqua_obj.pen.CapEx_pen / 1000000), '[Million $]')
    print(' '*2, "Feedbarge CapEx     ", "{:10.3f}".format(aqua_obj.pen.CapEx_feedbarge / 1000000), '[Million $]')

    print(' '*2, "fish feed OpEx      ", "{:10.3f}".format(aqua_obj.pen.fish_feed_price_annual / 1000000), '[Million $]')
    print(' '*2, "fingerling OpEx     ", "{:10.3f}".format(aqua_obj.pen.fingerling_price_annual / 1000000), '[Million $]')

    print(' '*2, "vessel travel NPV   ", "{:10.3f}".format(aqua_obj.vessel.cost_NPV / 1000000), '[Million $]', "{:2.1f}".format(100*aqua_obj.vessel.cost_NPV/aqua_obj.cost_NPV), '[%]')

    print(' '*2, "energy storage NPV  ", "{:10.3f}".format(aqua_obj.es.cost_NPV / 1000000), '[Million $]', "{:2.1f}".format(100*aqua_obj.es.cost_NPV/aqua_obj.cost_NPV), '[%]')

    print(' '*2, "diesel gen NPV      ", "{:10.3f}".format(aqua_obj.dieselgen.cost_NPV / 1000000), '[Million $]') #, "{:2.1f}".format(100*aqua_obj.dieselgen.cost_NPV/aqua_obj.cost_NPV), '[%]')
    print(' '*2, "diesel gen CapEx    ", "{:10.3f}".format(aqua_obj.dieselgen.CapEx / 1000000), '[Million $]')
    print(' '*2, "diesel OpEx         ", "{:10.3f}".format(aqua_obj.dieselgen.OpEx / 1000000), '[Million $]')
    print("-"*40)
    
def print_ineq_cons(title,aqua_obj):
    print_bold(title+" constraints:")
    print(' '*2, "normalized P_gen_cons              ", "{:10.3f}".format(aqua_obj.P_gen_cons), '[-]')
    print(' '*2, "normalized fish_yield_cons         ", "{:10.3f}".format(aqua_obj.fish_yield_cons), '[-]')
    print(' '*2, "normalized pen_ratio_low_cons      ", "{:10.3f}".format(aqua_obj.pen_ratio_low_cons), '[-]')
    print(' '*2, "normalized pen_ratio_up_cons       ", "{:10.3f}".format(aqua_obj.pen_ratio_up_cons), '[-]')
    print(' '*2, "sustainable_power_operation_cons   ", "{:10.3f}".format(aqua_obj.sustainable_power_operation_cons), '[-]')
    print("-"*40)

def print_location(title, aqua_obj):
    print_bold(title+" Location:")
    print(' '*2, "Longitude   ", "{:10.3f}".format(aqua_obj.pen.pos_long))
    print(' '*2, "Latitude    ", "{:10.3f}".format(aqua_obj.pen.pos_lat))
    print("-"*40)

def print_energy_storage(title, aqua_obj):
    print_bold(title+" es size:")
    print(' '*2, "es size   ", "{:10.3f}".format(aqua_obj.es.total_size), '[kWh]')
    print("-"*40)

def print_feedbarge(title, aqua_obj):
    print_bold(title+" feedbarge:")
    print(' '*2, "fish feed harvest week         ", "{:10.3f}".format(aqua_obj.pen.fish_feed_harvest_week), '[kg]')
    print(' '*2, "fish feed harvest week per pen ", "{:10.3f}".format(aqua_obj.pen.fish_feed_harvest_week / aqua_obj.pen.pen_number), '[kg]')
    print(' '*2, "feedbarge number               ", "{:10.3f}".format(aqua_obj.pen.feedbarge_number), '[-]')

def init_result(aqua_init_obj):
    print('+'*40) 
    if (aqua_init_obj.valid_point):
        print_objective("Initial",aqua_init_obj)
        print_P_rated("Initial",aqua_init_obj)
        print_price_breakdown("Initial",aqua_init_obj)
        print_ineq_cons("Initial",aqua_init_obj)
        print_carrying_capacity("Initial",aqua_init_obj)
    else:
        print('invalid init point')
    print('+'*40) 

def optimize_result(aqua_opt_obj, x_list, x_unit, res_opt):
    if (aqua_opt_obj.valid_point):
        print('optimization success: ',res_opt.success)
        col_width = len(max(x_list, key=len))
        
        print_bold("optimal design variable:")
        for i in range(len(x_list)):
            print(' '*2, x_list[i], ' '*(col_width - len(x_list[i])) , "{:10.3f}".format(res_opt.x[i]), x_unit[i])

        print(' '*2, "pen number    ", "{:10.3f}".format(aqua_opt_obj.pen.pen_number), '[-]')
        print("-"*40)
        
        print_location("deployment",aqua_opt_obj)
        print_objective("optimal",aqua_opt_obj)
        print_P_rated("optimal",aqua_opt_obj)
        print_price_breakdown("optimal",aqua_opt_obj)
        print_ineq_cons("optimal",aqua_opt_obj)
        print_carrying_capacity("optimal",aqua_opt_obj)
        print_feedbarge("optimal",aqua_opt_obj)
        print_energy_storage("optimal",aqua_opt_obj)
    else:
        print('invalid optimized point')
    print('+'*40) 

def bruteforce_result(title, x_best, res_opt, p):
    print('+'*80) 
    print_bold("brute force " + title)
    print('+'*80) 
    if x_best == None:
        print('invalid optimized point')
    else:
        if (res_opt.valid_point):
            col_width = len(max(x_best.list, key=len))
            
            print_bold("design variable:")
            for i in range(len(x_best.list)):
                print(' '*2, x_best.list[i], ' '*(col_width - len(x_best.list[i])) , "{:10.3f}".format(res_opt.x0[i]), x_best.unit[i])
            print("-"*40)
            
            print_objective("",res_opt)
            print_P_rated("",res_opt)
            print_price_breakdown("",res_opt)
        else:
            print('invalid optimized point')
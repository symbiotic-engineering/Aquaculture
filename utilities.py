import numpy as np
import pandas as pd
from modules import WPAF

def print_bold(str):
    print('\033[1m' + str + '\033[0;0m')
    return

def print_objective(title, wpaf):
    print_bold(title+" objective function terms:")
    print(' '*2, "Objective_func       ", "{:10.3f}".format(wpaf.obj_func))
    print(' '*2, "cost per yield       ", "{:10.3f}".format(wpaf.cost_per_yield), '[$ / kg]')
    print(' '*2, "NPV                  ", "{:10.3f}".format(wpaf.cost_NPV / 1000000), '[Million $]')
    print(' '*2, "levelized fish yield ", "{:10.3f}".format(wpaf.levelized_fish_yield / 1000000), '[kilo Tonne]')
    print(' '*2, "annual fish yield    ", "{:10.3f}".format(wpaf.aqua.fish_yield / 1000000), '[kilo Tonne]')
    print("-"*40)

def print_P_rated(title, wpaf):
    print_bold(title+" wave energy converter:")
    print(' '*2, "wec number    ", "{:10.3f}".format(wpaf.wec.wec_number), '[-]')
    #print(' '*2, "wave power    ", "{:10.3f}".format(wpaf.wave_in.P_wave/1000), '[kW/m]')
    print(' '*2, "wec_P_ave     ", "{:10.3f}".format(wpaf.wec.P_ave), '[kW]')
    print(' '*2, "wec AEP       ", "{:10.3f}".format(wpaf.wec.AEP), '[kWh]')
    print("-"*40)

def print_carrying_capacity(title, wpaf):
    print_bold(title+" carrying capacity:")
    print(' '*2, "TPF_O2              ", "{:10.3f}".format(wpaf.aqua.TPF_O2), '[kg fish/year]')
    print(' '*2, "Carrying Capacity   ", "{:10.3f}".format(wpaf.aqua.carrying_capacity), '[kg fish]')
    print("-"*40)

def print_price_breakdown(title, wpaf):
    print_bold(title+" price break down WPAF:")
    print(' '*2, "WPAF NPV                   ", "{:10.3f}".format(wpaf.cost_NPV / 1000000), '[Million $]',                              "{:5.1f}".format(100*wpaf.cost_NPV/wpaf.cost_NPV), '[%]')

    print(' '*2, "|__ wec NPV                ", "{:10.3f}".format(wpaf.wec.cost_NPV / 1000000), '[Million $]',                          "{:5.1f}".format(100*wpaf.wec.cost_NPV/wpaf.cost_NPV), '[%]')
    print(' '*2, "|   |__ wec CapEx          ", "{:10.3f}".format(wpaf.wec.CapEx / 1000000), '[Million $]',                             "{:5.1f}".format(100*wpaf.wec.CapEx/wpaf.cost_NPV), '[%]')
    print(' '*2, "|   |__ wec OpEx           ", "{:10.3f}".format(wpaf.wec.NPV_OpEx / 1000000), '[Million $]',                          "{:5.1f}".format(100*wpaf.wec.NPV_OpEx/wpaf.cost_NPV), '[%]')

    print(' '*2, "|__ aqua NPV               ", "{:10.3f}".format(wpaf.aqua.cost_NPV / 1000000), '[Million $]',                         "{:5.1f}".format(100*wpaf.aqua.cost_NPV/wpaf.cost_NPV), '[%]')
    print(' '*2, "|   |__ net pen CapEx      ", "{:10.3f}".format(wpaf.aqua.CapEx_netpen / 1000000), '[Million $]',                     "{:5.1f}".format(100*wpaf.aqua.CapEx_netpen/wpaf.cost_NPV), '[%]')
    print(' '*2, "|   |__ feedbarge CapEx    ", "{:10.3f}".format(wpaf.aqua.CapEx_feedbarge / 1000000), '[Million $]',                  "{:5.1f}".format(100*wpaf.aqua.CapEx_feedbarge/wpaf.cost_NPV), '[%]')

    print(' '*2, "|   |__ fish feed OpEx     ", "{:10.3f}".format(wpaf.aqua.NPV_OpEx_fish_feed_price / 1000000), '[Million $]',         "{:5.1f}".format(100*wpaf.aqua.NPV_OpEx_fish_feed_price/wpaf.cost_NPV), '[%]')
    print(' '*2, "|   |__ fingerling OpEx    ", "{:10.3f}".format(wpaf.aqua.NPV_OpEx_fingerling_price_annual / 1000000), '[Million $]', "{:5.1f}".format(100*wpaf.aqua.NPV_OpEx_fingerling_price_annual/wpaf.cost_NPV), '[%]')

    print(' '*2, "|   |__ vessel travel OpEx ", "{:10.3f}".format(wpaf.aqua.vessel.NPV_OpEx / 1000000), '[Million $]',                  "{:5.1f}".format(100*wpaf.aqua.vessel.NPV_OpEx/wpaf.cost_NPV), '[%]')

    print(' '*2, "|__ energy storage NPV     ", "{:10.3f}".format(wpaf.es.cost_NPV / 1000000), '[Million $]',                           "{:5.1f}".format(100*wpaf.es.cost_NPV/wpaf.cost_NPV), '[%]')
    print("-"*40)

    print_bold(title+" price break down DGAF:")
    print(' '*2, "DGAF NPV                   ", "{:10.3f}".format(wpaf.cost_NPV_diesel / 1000000), '[Million $]',                         "{:5.1f}".format(100*wpaf.cost_NPV_diesel/wpaf.cost_NPV_diesel), '[%]')

    print(' '*2, "|__ aqua NPV               ", "{:10.3f}".format(wpaf.aqua.cost_NPV / 1000000), '[Million $]',                         "{:5.1f}".format(100*wpaf.aqua.cost_NPV/wpaf.cost_NPV_diesel), '[%]')
    print(' '*2, "|   |__ net pen CapEx      ", "{:10.3f}".format(wpaf.aqua.CapEx_netpen / 1000000), '[Million $]',                     "{:5.1f}".format(100*wpaf.aqua.CapEx_netpen/wpaf.cost_NPV_diesel), '[%]')
    print(' '*2, "|   |__ feedbarge CapEx    ", "{:10.3f}".format(wpaf.aqua.CapEx_feedbarge / 1000000), '[Million $]',                  "{:5.1f}".format(100*wpaf.aqua.CapEx_feedbarge/wpaf.cost_NPV_diesel), '[%]')

    print(' '*2, "|   |__ fish feed OpEx     ", "{:10.3f}".format(wpaf.aqua.NPV_OpEx_fish_feed_price / 1000000), '[Million $]',         "{:5.1f}".format(100*wpaf.aqua.NPV_OpEx_fish_feed_price/wpaf.cost_NPV_diesel), '[%]')
    print(' '*2, "|   |__ fingerling OpEx    ", "{:10.3f}".format(wpaf.aqua.NPV_OpEx_fingerling_price_annual / 1000000), '[Million $]', "{:5.1f}".format(100*wpaf.aqua.NPV_OpEx_fingerling_price_annual/wpaf.cost_NPV_diesel), '[%]')

    print(' '*2, "|   |__ vessel travel OpEx ", "{:10.3f}".format(wpaf.aqua.vessel.NPV_OpEx / 1000000), '[Million $]',                  "{:5.1f}".format(100*wpaf.aqua.vessel.NPV_OpEx/wpaf.cost_NPV_diesel), '[%]')

    print(' '*2, "|__ diesel gen NPV         ", "{:10.3f}".format(wpaf.dieselgen.cost_NPV / 1000000), '[Million $]',                    "{:5.1f}".format(100*wpaf.dieselgen.cost_NPV/wpaf.cost_NPV_diesel), '[%]')
    print(' '*2, "    |__ diesel gen CapEx   ", "{:10.3f}".format(wpaf.dieselgen.CapEx / 1000000), '[Million $]',                       "{:5.1f}".format(100*wpaf.dieselgen.CapEx/wpaf.cost_NPV_diesel), '[%]')
    print(' '*2, "    |__ diesel OpEx        ", "{:10.3f}".format(wpaf.dieselgen.NPV_OpEx / 1000000), '[Million $]',                    "{:5.1f}".format(100*wpaf.dieselgen.NPV_OpEx/wpaf.cost_NPV_diesel), '[%]')
    print("-"*40)
    
def print_ineq_cons(title,wpaf):
    print_bold(title+" constraints:")
    #print(' '*2, "normalized P_gen_cons              ", "{:10.3f}".format(wpaf.P_gen_cons), '[-]')
    print(' '*2, "normalized fish_yield_cons         ", "{:10.3f}".format(wpaf.fish_yield_cons), '[-]')
    # print(' '*2, "normalized pen_ratio_low_cons      ", "{:10.3f}".format(wpaf.pen_ratio_low_cons), '[-]')
    # print(' '*2, "normalized pen_ratio_up_cons       ", "{:10.3f}".format(wpaf.pen_ratio_up_cons), '[-]')
    print(' '*2, "sustainable_power_operation_cons   ", "{:10.3f}".format(wpaf.sustainable_power_operation_cons), '[-]')
    print("-"*40)

def print_location(title, wpaf):
    print_bold(title+" Deployment Site:")
    print(' '*2, "Longitude               ", "{:10.3f}".format(wpaf.aqua.pos_long))
    print(' '*2, "Latitude                ", "{:10.3f}".format(wpaf.aqua.pos_lat))
    print(' '*2, "Hs annual avg           ", "{:10.3f}".format(wpaf.wec.wave.Hs.mean()), '[m]')
    print(' '*2, "Tp annual avg           ", "{:10.3f}".format(wpaf.wec.wave.Tp.mean()), '[s]')
    print(' '*2, "Temp annual avg         ", "{:10.3f}".format(wpaf.aqua.temp.mean()), '[C]')
    print(' '*2, "Bathymetry              ", "{:10.3f}".format(-wpaf.aqua.bathymetry), '[m]')
    print(' '*2, "DO2 annual avg          ", "{:10.3f}".format(wpaf.aqua.O2_in), '[mg/l]')
    print(' '*2, "Salinity annual avg     ", "{:10.3f}".format(wpaf.aqua.salinity), '[PSU]')
    print(' '*2, "Current Speed annual avg", "{:10.3f}".format(wpaf.aqua.U), '[m/s]')
    print("-"*40)

def print_energy_storage(title, wpaf):
    print_bold(title+" es size:")
    print(' '*2, "es size   ", "{:10.3f}".format(wpaf.es.total_size), '[kWh]')
    print("-"*40)

def print_netpen_parameter(title, wpaf):
    print_bold(title+" netpen:")
    print(' '*2, "pen number    ", "{:10.3f}".format(wpaf.aqua.pen_number), '[-]')
    print(' '*2, "pen height    ", "{:10.3f}".format(wpaf.aqua.netpen.H), '[kWh]')
    print("-"*40)

def print_feedbarge(title, wpaf):
    print_bold(title+" feedbarge:")
    print(' '*2, "fish feed harvest week         ", "{:10.3f}".format(wpaf.aqua.fish_feed_harvest_week), '[kg]')
    print(' '*2, "fish feed harvest week per pen ", "{:10.3f}".format(wpaf.aqua.fish_feed_harvest_week / wpaf.aqua.pen_number), '[kg]')
    print(' '*2, "feedbarge number               ", "{:10.3f}".format(wpaf.aqua.feedbarge_number), '[-]')

def print_init_result(aqua_init_obj):
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

def print_soo_optimize_result(wpaf_opt_obj, x_list, x_unit, res_opt):
    if (wpaf_opt_obj.valid_point):
        print('optimization success: ',res_opt.success)
        col_width = len(max(x_list, key=len))
        
        print_bold("optimal design variable:")
        for i in range(len(x_list)):
            print(' '*2, x_list[i], ' '*(col_width - len(x_list[i])) , "{:10.3f}".format(res_opt.x[i]), x_unit[i])

        
        print("-"*40)
        
        print_objective("optimal",wpaf_opt_obj)
        print_P_rated("optimal",wpaf_opt_obj)
        print_price_breakdown("optimal",wpaf_opt_obj)
        print_ineq_cons("optimal",wpaf_opt_obj)
        print_netpen_parameter("optimal",wpaf_opt_obj)
        print_location("deployment",wpaf_opt_obj)
        print_carrying_capacity("optimal",wpaf_opt_obj)
        print_feedbarge("optimal",wpaf_opt_obj)
        print_energy_storage("optimal",wpaf_opt_obj)
    else:
        print('invalid optimized point')
    print('+'*40) 

def print_moo_optimize_result(wpaf_opt_obj, x_list, x_unit, moo_res_opt):
    if (wpaf_opt_obj.valid_point):
        if moo_res_opt is None:
            print('optimization success: False')
        else:
            print('optimization success: True')

        col_width = len(max(x_list, key=len))
        
        print_bold("optimal design variable:")
        for i in range(len(x_list)):
            try:
                print(' '*2, x_list[i], ' '*(col_width - len(x_list[i])) , "{:10.3f}".format(moo_res_opt.X[i]), x_unit[i])
            except:
                print(' '*2, x_list[i], ' '*(col_width - len(x_list[i])) , "{:10.3f}".format(moo_res_opt.X[0][i]), x_unit[i])

        print("-"*40)
        
        print_objective("optimal",wpaf_opt_obj)
        print_P_rated("optimal",wpaf_opt_obj)
        print_price_breakdown("optimal",wpaf_opt_obj)
        print_ineq_cons("optimal",wpaf_opt_obj)
        print_netpen_parameter("optimal",wpaf_opt_obj)
        print_location("deployment",wpaf_opt_obj)
        print_carrying_capacity("optimal",wpaf_opt_obj)
        print_feedbarge("optimal",wpaf_opt_obj)
        print_energy_storage("optimal",wpaf_opt_obj)
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
import numpy as np
import copy
from scipy.optimize import minimize
from functools import partial
import modules

class OpData():
    def __init__(self, name):
        self.name = name
        self.list, self.nom_dict, self.unit, self.bnds, self.label = default_value(name)
    
    @property
    def nom0(self):
        data0 = []
        for i in range(len(self.list)):
            data0.append(self.nom_dict[self.list[i]])
        return data0


class OpObj(object):
    def __init__(self, x0, x_name, p, max_iter):
        self.x_name, self.p = x_name, p
        self.x0 = x0
        self.f = np.full(shape=(max_iter+1,), fill_value=np.NaN)
        self.ineq = np.full(shape=(max_iter+1,len(modules.ineq_constraint(x0, x_name, p))), fill_value=np.NaN)
        self.eq = np.full(shape=(max_iter+1,0), fill_value=np.NaN)
        self.count = 0
        
    def obj_fun(self, x):
        return obj_fun(x, self.x_name, self.p)

def cb(xk, obj=None):
    obj.f[obj.count] = obj.obj_fun(xk)
    obj.ineq[obj.count] = modules.ineq_constraint(xk, obj.x_name, obj.p)
    obj.eq[obj.count] = modules.eq_constraint(xk, obj.x_name, obj.p)
    obj.count += 1

def obj_fun(x0, x_name, p):
    return modules.obj(x0, x_name, p)

def default_value(v_name):
    v_label = ''
    for i in range(len(v_name)):
        if (i!=0):
            v_label += ' & '
        v_label += v_name[i]
    v_list = modules.variable_lookup(v_name)
    v_list_default_values = modules.default_values(v_name)
    v_list_bnds_values = modules.bnds_values(v_name)
    v_nom = {}
    v_unit = []
    v_bnds = []
    for i in range(len(v_list)):
        v_nom[v_list[i]] = v_list_default_values[v_list[i]][0]
        v_unit.append(v_list_default_values[v_list[i]][1])
        if v_list[i] in v_list_bnds_values.keys():
            v_bnds.append(v_list_bnds_values[v_list[i]])
    return v_list, v_nom, v_unit, v_bnds, v_label


def argument_fun(x_name, p_name, p_vals, all_vars):
    all_input = x_name + p_name
    default_vars = []
    for i in range(len(all_vars)):
        if all_vars[i] not in all_input:
            default_vars.append(all_vars[i])
    p = OpData(default_vars)
    
    # fill non-default parameters
    if p_name!=[]:
        p.name = p.name + p_name
        new_list = modules.variable_lookup(p_name)
        p.list = p.list + new_list
        for i in range(len(new_list)):
            p.nom_dict[new_list[i]] = p_vals[new_list[i]]
          
    return p

#'Finding optimal '+ x + ' while holding '+ p + ' constant.'
def run_optimization(x_name, x_vals, p_name, p_vals, all_vars, max_iter):
    # optimizes the design variables x_name, with parameters p_name set to 
    # non-default values p_vals, and other parameters set to default values.
    
    # design variables
    x = OpData(x_name)
    #x = [x_data.label,x_data.list]
    x0 = []
    if x_vals==[]:
        x0 = x.num0
    else:
        x0 = x_vals
    
    # fill default parameters
    p = argument_fun(x.name, p_name, p_vals, all_vars)

        
    # set up optimization problem
    op_obj = OpObj(x0, x.name, p.nom_dict, max_iter) 
    arguments = (x.name, p.nom_dict)
    cons = []
    cons.append({'type': 'ineq', 
                 'fun': modules.ineq_constraint, 'args': arguments})
    
    # adding x_bounds as constrains for COBYLA
    for factor in range(len(x.bnds)):
        lower, upper = x.bnds[factor]
        l = {'type': 'ineq',
             'fun': lambda x, lb=lower, i=factor: x[i] - lb}
        u = {'type': 'ineq',
             'fun': lambda x, ub=upper, i=factor: ub - x[i]}
        cons.append(l)
        cons.append(u)

    options={"maxiter":max_iter}
    
    res = minimize(obj_fun, op_obj.x0, 
                   args=arguments, 
                   #method='COBYLA',
                   method='SLSQP',
                   #bounds=x_bnds, 
                   constraints=cons,
                   options=options,
                   callback=partial(cb, obj=op_obj))
    
    return res, op_obj, p
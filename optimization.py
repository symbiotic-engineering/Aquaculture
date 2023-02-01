import numpy as np
from scipy.optimize import minimize
from functools import partial
import modules


class OpObj(object):
    def __init__(self, x0, x_name, p, max_iter):
        self.x_name, self.p = x_name, p
        self.x0 = x0
        self.f = np.full(shape=(max_iter,), fill_value=np.NaN)
        self.x_history = np.full(shape=(max_iter,len(x0)), fill_value=np.NaN)
        self.obj_history = np.full(shape=(max_iter,len(modules.obj_terms(x0, x_name, p))), fill_value=np.NaN)
        self.ineq = np.full(shape=(max_iter,len(modules.ineq_constraint(x0, x_name, p))), fill_value=np.NaN)
        self.eq = np.full(shape=(max_iter,0), fill_value=np.NaN)
        self.count = 0
        
    def obj_fun(self, x):
        return obj_fun(x, self.x_name, self.p)

def cb(xk, obj=None):
    obj.f[obj.count] = obj.obj_fun(xk)
    obj.x_history[obj.count] = xk
    obj.obj_history[obj.count] = modules.obj_terms(xk, obj.x_name, obj.p)
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
    p_list, p, p_unit, p_bnds, p_label = default_value(default_vars)
    
    # fill non-default parameters
    if p_name!=[]:
        p_list = modules.variable_lookup(p_name)
        for i in range(len(p_list)):
            p[p_list[i]] = p_vals[p_list[i]]
          
    return p_list, p, p_unit, p_bnds, p_label

#'Finding optimal '+ x + ' while holding '+ p + ' constant.'
def run_optimization(x_name, x_vals, p_name, p_vals, all_vars, max_iter):
    # optimizes the design variables x_name, with parameters p_name set to 
    # non-default values p_vals, and other parameters set to default values.
    
    # design variables
    x_list, x0_val, x_unit, x_bnds, x_label = default_value(x_name)
    x = [x_label,x_list]
    x0 = []
    if x_vals==[]:
        for i in range(len(x_list)):
            x0.append(x0_val[x_list[i]])
    else:
        x0 = x_vals
    
    # fill default parameters
    p_list, p, p_unit, p_bnds, p_label = argument_fun(x_name, p_name, p_vals, all_vars)

        
    # set up optimization problem
    op_obj = OpObj(x0, x_name, p, max_iter) 
    arguments = (x_name, p)
    cons = []
    cons.append({'type': 'ineq', 
                 'fun': modules.ineq_constraint, 'args': arguments})
    
    # adding x_bounds as constrains for COBYLA
    for factor in range(len(x_bnds)):
        lower, upper = x_bnds[factor]
        l = {'type': 'ineq',
             'fun': lambda x, lb=lower, i=factor: x[i] - lb}
        u = {'type': 'ineq',
             'fun': lambda x, ub=upper, i=factor: ub - x[i]}
        cons.append(l)
        cons.append(u)

    options={"maxiter":max_iter}  #"ftol": 1e-4
    
    
    res = minimize(obj_fun, op_obj.x0, 
                   args=arguments, 
                   method='SLSQP',
                   #bounds=x_bnds, 
                   constraints=cons,
                   options=options,
                   callback=partial(cb, obj=op_obj))
    
    return res, op_obj
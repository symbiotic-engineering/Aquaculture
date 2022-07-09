ratios = [.8 1 1.2];
% sensitivity of optimal wec and pen design to environment/location
dx_wec_pen_star_dx_env = sensitivity_sweep('x_env', ratios, {'x_wec','x_pen'});

% sensitivity of optimal environment/location to wec and pen design
dx_env_dx_wec_pen = sensitivity_sweep({'x_wec','x_pen'}, ratios, 'x_env');

function dv2_star_dv1 = sensitivity_sweep(v1_name, v1_ratios, v2_name)
% v1: variable being swept, held constant for each optimization
% v2: variable being optimized

    v1_list = variable_lookup(v1_name);
    v2_list = variable_lookup(v2_name);
    v1_nom = default_values(v1_name);
    num_v1 = length(v1_list);
    num_v2 = length(v2_list);
    num_ratios = length(v1_ratios);
    
    v2_star = zeros(num_v1,num_ratios,num_v2);
    dv2_star_dv1 = zeros(num_v1,num_v2);
    
    for i = 1:num_v1
        for j = 1:num_ratios
            v1 = v1_nom;
            v1.(v1_list{i}) = v1_nom.(v1_list{i}) * v1_ratios(j);
            v2_star(i,j,:) = run_optimization(v2_name,v1_name,v1);
        end
        idx_nom = ratios==1;
        dv2_star = (v2_star(i,end,:) - v2_star(i,1,:)) ./ v2_star(i,idx_nom,:);
        dv2_star_dv1(i,:) = dv2_star / (v1_ratios(end) - v1_ratios(1));
    end

end

function x_star = run_optimization(x_name,p_name,p_vals)
% optimizes the design variables x_name, with parameters p_name set to 
% non-default values p_vals, and other parameters set to default values.

    if nargin==1
        p_name = {};
    end

    % create optimization variable
    x_list = variable_lookup(x_name);
    x_label = [x_name{:}];
    x = optimvar(x_label,x_list);

    % create struct to hold optimization variable
    for i = 1:length(x_list)
        x_struct.(x_list{i}) = x(x_list{i});
    end
    
    % fill default parameters
    all_vars = {'x_wec','x_pen','x_env','x_fish','x_consts'};
    is_default = ~contains(all_vars,[x_name p_name]);
    default_vars = all_vars(is_default);
    p = default_values(default_vars);

    % fill non-default parameters
    if nargin>1
        p_list = variable_lookup(p_name);
        for i = 1:length(p_list)
            p.(p_list{i}) = p_vals.(p_list{i});
        end
    end

    % set up optimization problem
    [J,g,h] = simulate(x_struct,p);
    desc = ['Finding optimal ' x_name ' while holding ' p_name ' constant.'];
    problem = optimproblem('Objective',J,'Description',desc);
    problem.constraints.ineq = g >= 0;
    problem.constraints.eq   = h == 0;
    
    % solve optimization problem
    x_star = solve(problem);

end

function vals = default_values(var_category_names)
    if any(strcmp('x_wec',var_category_names))
        vals.wec_type = 1;
        vals.capture_width = 1;
    end
    if any(strcmp('x_pen',var_category_names))
        vals.pen_diameter = 1;
        vals.pen_height = 1;
        vals.num_pens = 1;
        vals.spacing = 1;
        vals.stocking_density = 1;
    end
    if any(strcmp('x_env',var_category_names))
        vals.temp = 1;
        vals.salinity = 1;
        vals.O2_in = 1;
        vals.min_current_speed = 1;
        vals.wave_height = 1;
        vals.wave_period = 1;
    end
    if any(strcmp('x_consts',var_category_names))
        vals.wec_unit_cost = 1;
        vals.pen_unit_cost = 1;
        vals.permeability = 1;
        vals.capture_width_ratio_dict = 1;
        vals.wave_damping_dict = 1;
    end
    if any(strcmp('x_fish',var_category_names))
        vals.F_p = 1;
        vals.F_c = 1;
        vals.F_f = 1;
        vals.A_p = 1;
        vals.A_c = 1;
        vals.A_f = 1;
        vals.O_p = 1;
        vals.O_c = 1;
        vals.O_f = 1;
        vals.C_p = 1;
        vals.C_c = 1;
        vals.C_f = 1;
        vals.P_f = 1;
        vals.P_p = 1;
        vals.O2_min = 1;
        vals.tau = 1;
        vals.loss_rate = 1;
        vals.harvest_weight = 1;
    end

    %assert(fieldnames(vals) == variable_lookup(var_category_names));
end

function var_list = variable_lookup(var_category_names)
% input is a cell array containing some subset of the following categories:
% x_wec, x_pen, x_env, x_fish
% output is a cell array with variable names matching the input categoeies.

    var_list = {};
    if any(strcmp('x_wec',var_category_names))
        var_list = [var_list, {'wec_type','capture_width'}];
    end
    if any(strcmp('x_pen',var_category_names))
        var_list = [var_list, {'pen_diameter','pen_height','num_pens',...
                               'spacing','stocking_density'}];
    end
    if any(strcmp('x_env',var_category_names))
        var_list = [var_list, {'temp','O2_in','min_current_speed',...
                               'wave_height','wave_period'}];
    end
    if any(strcmp('x_consts',var_category_names))
        var_list = [var_list, {'wec_unit_cost','pen_unit_cost','permeability',...
                            'capture_width_ratio_dict','wave_damping_dict'}];
    end
    if any(strcmp('x_fish',var_category_names))
        var_list = [var_list, {'F_f','F_p','F_c',...
                               'A_f','A_p','A_c',...
                               'O_f','O_p','O_c',...
                               'C_f','C_p','C_c','P_f','P_p',...
                               'O2_min','tau','loss_rate','harvest_weight'}];
    end
    if isempty(var_list)
        error('Your input did not match any of the category names.')
    end

end

function [J,g,h] = simulate(x,p)
%SIM Summary of this function goes here
%   Detailed explanation goes here
    mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
    ins = mergestructs(x,p);
    
    wec = WEC(ins.capture_width, ins.capture_width_ratio_dict, ...
            ins.wave_damping_dict, ins.wec_type, ins.wec_unit_cost);

    wave_in = Wave(ins.wave_height, ins.wave_period);

    pen = Pen(ins.pen_diameter, ins.pen_height, ins.stocking_density, ...
            ins.num_pens, ins.spacing, ins.pen_unit_cost, ins.loss_rate, ...
            ins.harvest_weight, ins.env_params);

    % run each module 
    wave_out = wave_climate(wec,wave_in);
    fish_yield = fish(wave_out,pen);
    pow = power(wec, wave_in);
    carrying_capacity = environment(pen);
    price = econ(wec, pen);

    % outputs
    cost_per_yield = price/fish_yield;
    J = cost_per_yield;
    g = [pow, carrying_capacity];
    h = [];

end

function P_gen = power(wec, wave)
    %assert(isinstance(wec,WEC))
    %assert(isinstance(wave,Wave))

    P_gen = wave.power * wec.capture_width * wec.capture_width_ratio;
end

function price = econ(wec, pen)
    %assert(isinstance(wec,WEC))
    %assert(isinstance(pen,Pen))

    price = wec.price + pen.price;
end

function wave_out = wave_climate(wec, wave)
    %assert(isinstance(wec,WEC))
    %assert(isinstance(wave,Wave))

    Hs = wec.wave_damping * wave.Hs;
    T = wave.T;
    wave_out = Wave(Hs,T);
end

function fish_yield = fish(wave, pen)
    %assert(isinstance(wave,Wave))
    %assert(isinstance(pen,Pen))

    if wave.Hs > 3
        extra_loss_rate = 0.1;
    else
        extra_loss_rate = 0;
    end

    survival_rate = (1-(pen.loss_rate + extra_loss_rate));
    fish_yield = pen.SD * survival_rate * pen.volume * pen.harvest_weight;
end

function cc = environment(pen)
    %assert(isinstance(pen,Pen))
    cc = pen.carrying_capacity;
end

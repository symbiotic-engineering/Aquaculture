classdef Pen
    %PENCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        D
        H
        SD
        n
        spacing
        unit_cost
        loss_rate
        harvest_weight

        O2_in
        O2_min
        P_f
        P_p
        U_min
        tau
        permeability

        F_f
        F_p
        F_c

        A_f
        A_p
        A_c

        O_f
        O_p
        O_c

        C_f
        C_p
        C_c
    end
    properties (Dependent)
        DO2
        carrying_capacity
        price
        volume
    end
    
    methods
        function obj = Pen(D,H,SD,n,spacing,unit_cost,loss_rate,harvest_weight,env_params)
            obj.D = D; 
            obj.H = H;
            obj.SD = SD;
            obj.n = n;
            obj.spacing = spacing;
    
            obj.unit_cost = unit_cost;
            obj.loss_rate = loss_rate;
            obj.harvest_weight = harvest_weight;
    
            obj.O2_in = env_params.O2_in;
            obj.O2_min = env_params.O2_min;
            obj.P_f = env_params.P_f;
            obj.P_p = env_params.P_p;
            obj.U_min = env_params.U_min;
            obj.tau = env_params.tau;
            obj.permeability = env_params.permeability;
    
            obj.F_f = env_params.F_f;
            obj.F_p = env_params.F_p;
            obj.F_c = env_params.F_c;
    
            obj.A_f = env_params.A_f;
            obj.A_p = env_params.A_p;
            obj.A_c = env_params.A_c;
            
            obj.O_f = env_params.O_f;
            obj.O_p = env_params.O_p;
            obj.O_c = env_params.O_c;
            
            obj.C_f = env_params.C_f;
            obj.C_p = env_params.C_p;
            obj.C_c = env_params.C_c;
        end
        
        function DO2 = get.DO2(obj)
            % specific energy content of feed
            delta = obj.F_f * obj.C_f + obj.F_p * obj.C_p + obj.F_c * obj.C_c;
    
            % specific energy content of fish
            C_f_star = 0.85 * obj.C_p * obj.P_p + obj.C_f * obj.P_f;
    
            % fraction of food energy contributed by progeins, fat, and carbs
            E_p = obj.F_p * obj.C_p / delta;
            E_f = obj.F_f * obj.C_f / delta;
            E_c = obj.F_c * obj.C_c / delta;
    
            % metabolizable energy content of food
            FL = (1-obj.A_p) * E_p + (1-obj.A_f) * E_f + (1-obj.A_c) * E_c;
            BC = 0.3 * obj.A_p * E_p + 0.05 * (obj.A_f * E_f + obj.A_c * E_c);
            eps = 1 - FL - BC;
            eps_star = eps - 0.15 * obj.F_p * obj.C_p * obj.A_p / delta;
    
            % Water temperature as a function of time
            time = linspace(0,51,1); % time vector [weeks]
            T = 52;                  % period [weeks]
            w = 2*pi/T;              % frequency [1/weeks]
            phi = 2*pi/3;            % phase offset [-]
            T_max = 23;
            T_min = 4;
            T_bar = (T_max+T_min)/2;
            T_amp = (T_max-T_min)/2;
            Temp = T_bar + T_amp * cos(w * time + phi);
    
            % Fish growth as a function of time
            a = 0.038;
            W_0 = 0;
            %integral = trapz( exp(Temp*obj.tau), x=time )
            integral = Temp; % fixme
            W = (W_0^(1/3) + a/3 * integral)^3;
    
            % Growth rate as a function of time
            b = 2/3;
            W_dot = a * W^b * exp(Temp * obj.tau);
    
            % Rate of energy ingested by fish, cal/day
            alpha = 11;
            gamma = 0.8;
            Q_r = 1/eps_star * (alpha * W^gamma + a * C_f_star * W^b) * exp(time*obj.tau);
    
            % Respiratory oxygen demand with respect to protein, fat, and carb consumption of fish
            DO2_p = (obj.F_p * obj.A_p * Q_r / delta - obj.P_p * W_dot) * obj.O_p;
            DO2_f = (obj.F_f * obj.A_f * Q_r / delta - obj.P_f * W_dot) * obj.O_f;
            DO2_c = obj.F_c * obj.A_c * Q_r / delta * obj.O_c;
    
            % Total respiratory oxygen demand of fish per day
            DO2 = DO2_p + DO2_f + DO2_c;
        end
        function cc = get.carrying_capacity(obj)
            length = obj.n * obj.D + obj.spacing * (obj.n-1);
            cc = (obj.O2_in - obj.O2_min) * length * obj.H * obj.permeability * obj.U_min / obj.DO2;
            cc = min(cc);
        end
        function price = get.price(obj)
            price = obj.D * obj.H * obj.unit_cost;
        end
        function volume = get.volume(obj)
            volume = pi * obj.D^2 / 4 * obj.H;
        end
    end
end


classdef WEC
    %WECCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        capture_width
        capture_width_ratio_dict
        wave_damping_dict
        wec_type
        unit_cost
    end
    properties (Dependent)
        price
        wave_damping
        capture_width_ratio
    end
    
    methods
        function obj = WEC(capture_width,capture_width_ratio_dict,wave_damping_dict,wec_type,unit_cost)
            %WECCLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.capture_width = capture_width;
            obj.capture_width_ratio_dict = capture_width_ratio_dict;
            obj.wave_damping_dict = wave_damping_dict;
            obj.wec_type = wec_type;
            obj.unit_cost = unit_cost;
        end
        
        function price = get.price(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            price = obj.capture_width * obj.unit_cost;
        end
        function wave_damping = get.wave_damping(obj)
            wave_damping = obj.capture_width_ratio_dict.(obj.wec_type);
        end
        function CWR = get.capture_width_ratio(obj)
            CWR = obj.capture_width_ratio_dict.(obj.wec_type);
        end
    end
end


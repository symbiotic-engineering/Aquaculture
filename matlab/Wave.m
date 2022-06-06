classdef Wave
    %WAVECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Hs
        T
        rho = 1000
        g = 9.81
    end
    
    methods
        function obj = Wave(Hs,T)
            %WAVECLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Hs = Hs;
            obj.T = T;
        end
        
        function power = power(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            power = 1/32 * 1/pi * obj.rho * obj.g^2 * obj.Hs^2 * obj.T;
        end
    end
end


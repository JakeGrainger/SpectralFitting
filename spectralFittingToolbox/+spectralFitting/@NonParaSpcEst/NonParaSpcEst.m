classdef NonParaSpcEst
    properties
        % computed spectral estimate
        I
        % angular frequencies corresponding to I
        omega
    end
    methods
        function obj = NonParaSpcEst(I,omega)
            obj.I = I;
            obj.omega = omega;
        end
        plot(obj,varargin)
        plotdB(obj,varargin)
    end
end
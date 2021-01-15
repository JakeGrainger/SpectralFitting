classdef Periodogram < spectralFitting.NonParaSpcEst
    methods
        function obj = Periodogram(timeseries,Delta)
            arguments
                timeseries (:,1) {mustBeNumeric}
                Delta(1,1) {mustBeNumeric}
            end
            N = length(timeseries);
            I = Delta/(N*2*pi)*abs(fft(timeseries)).^2;
            I = I(1:floor(end/2)+1);
            omega = 2*pi/Delta*(0:1/N:1/2)';
            obj@spectralFitting.NonParaSpcEst(I,omega);
        end
    end
end
    

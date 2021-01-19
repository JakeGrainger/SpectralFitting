classdef BartlettPeriodogram < spectralFitting.NonParaSpcEst
    methods
        function obj = BartlettPeriodogram(timeseries,Delta,windowSize)
            arguments
                timeseries (:,1) {mustBeNumeric}
                Delta(1,1) {mustBeNumeric}
                windowSize(1,1) {mustBeNumeric}
            end
            nSegments = floor(length(timeseries)/windowSize);
            segmentPeriodograms = nan(windowSize,nSegments);
            for iSeg = 1:nSegments
                segmentPeriodograms(:,iSeg) = Delta/(windowSize*2*pi)*...
                    abs(fft(timeseries((iSeg-1)*windowSize+1:iSeg*windowSize,:),[],1)).^2;
            end
            S = mean(segmentPeriodograms,2);
            I = S(1:floor(end/2)+1);
            omega = 2*pi/Delta*(0:1/windowSize:1/2)';
            obj@spectralFitting.NonParaSpcEst(I,omega);
        end
    end
end
    
 

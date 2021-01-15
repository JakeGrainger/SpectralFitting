function corrI = computeCorrPeriodogram(obj,N,Delta)
% computeCorrPeriodogram  Compute the correlation of the periodogram.
%   covI = computeCorrPeriodogram(obj,N,Delta) computes the correlation of
%   the periodogram at the Fourier frequencies defined by N and Delta.
%
%   See also computeCovPeriodogram, plotCovI, plotCorrI.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
end
covI = parCovPeriodogram(obj,N,Delta,obj.parameter);
corrI = corrcov(covI);
end

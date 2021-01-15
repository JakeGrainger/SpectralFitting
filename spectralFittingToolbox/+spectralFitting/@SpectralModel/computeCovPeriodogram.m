function covI = computeCovPeriodogram(obj,N,Delta)
% computeCovPeriodogram  Compute the covariance of the periodogram.
%   covI = computeCovPeriodogram(obj,N,Delta) computes the covariance of
%   the periodogram at the Fourier frequencies defined by N and Delta.
%
%   See also computeCorrPeriodogram, plotCovI, plotCorrI.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
end
covI = parCovPeriodogram(obj,N,Delta,obj.parameter);
end

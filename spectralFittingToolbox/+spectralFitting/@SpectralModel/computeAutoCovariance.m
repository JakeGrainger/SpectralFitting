function varargout = computeAutoCovariance(obj,N,Delta)
% computeAutoCovariance  Compute the autocovariance.
%   [C,dC,d2C] = computeAutoCovariance(obj,N,Delta) provides the
%   autocovariance and first and second derivatives, at lags given by
%   0,Delta,...,N-1*Delta. The second derivatives are in the same format as
%   computeSpectralDensity. If the property autoCovariance is left empty, 
%   then this will be approximated from the
%   aliased spectral density.
%
%   See also computeSpectralDensity, computeAliasedSpectralDensity,
%   computeExpectedPeriodogram, plotAcv.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
end

[varargout{1:max(1,nargout)}] =...
    parAutoCovariance(obj,N,Delta,obj.parameter);
end
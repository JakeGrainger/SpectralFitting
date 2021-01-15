function varargout = computeExpectedPeriodogram(obj,N,Delta)
% computeExpectedPeriodogram  Compute the expected periodogram.
%   [E,dE,d2E] = computeExpectedPeriodogram(obj,N,Delta) computes the
%   expected peridogram at Fourier frequencies specified by N and Delta. 
%   Derivatives and second derivatives given in the same format as
%   computeSpectralDensity.
%
%   See also computeSpectralDensity, computeAliasedSpectralDensity,
%   computeAutoCovariance, plotEI.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
end
[varargout{1:max(1,nargout)}] = parExpectedPeriodogram(obj,N,Delta,obj.parameter,1:floor(N/2)+2);
end
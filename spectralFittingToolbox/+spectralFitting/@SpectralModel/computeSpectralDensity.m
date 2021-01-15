function varargout = computeSpectralDensity(obj,omega,Delta)
% computeSpectralDensity  Compute the spectral density.
%   S = computeSpectralDensity(obj,omega,Delta) computes the spectral
%   density of the model defined in obj at the angular frequencies omega.
%   If differencing is enables, computes the differenced spectral density.
%   
%   [S,dS,d2S] = computeSpectralDensity(obj,omega,Delta) also computes
%   first and second derivatives. Note the second derivatives are only
%   computed for the parameter pairs specified to be non-zero by the
%   hessIndex property. By default, this is the lower triangle. The second
%   derivative will be a series of columns where each column is the
%   derivatives with respect to two parameters for each frequency. This
%   should be specified by the spectralDensity property set by the user.
%
%   See also computeAliasedSpectralDensity, computeAutoCovariance,
%   computeExpectedPeriodogram, plotSdf.
arguments
    obj spectralFitting.SpectralModel
    omega {mustBeNumeric,mustBeNonempty}
    Delta {mustBeNumeric} = [] % Delta only needs to be specified if differencing is being used.
end
if obj.useDifferencing
    if isempty(Delta)
        error("Please specify Delta when using differencing.")
    end
    if numel(Delta) > 1
        error("Delta must be scalar.")
    end
end
[varargout{1:max(1,nargout)}] = parSpectralDensity(obj,omega,Delta,obj.parameter);
end
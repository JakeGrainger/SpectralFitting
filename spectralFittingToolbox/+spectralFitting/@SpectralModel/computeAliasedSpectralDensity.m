function varargout = computeAliasedSpectralDensity(obj,omega,Delta)
% computeAliasedSpectralDensity  Compute the aliased spectral density.
%   S = computeAliasedSpectralDensity(obj,omega,Delta) computes the
%   aliased spectral density at the provided angular frequencies, omega. 
%   If differencing is enabled, will the aliased spectral density of the
%   differenced process. These will be approximated if the property 
%   aliasedSpectralDensity is empty.
%
%   [S,dS,d2S] = computeAliasedSpectralDensity(obj,omega,Delta) will
%   compute the derivatives as well. The second derivatives are in the same
%   format as computeSpectralDensity.
%
%   See also computeSpectralDensity, computeAutoCovariance,
%   computeExpectedPeriodogram, plotAliasedSdf.
arguments
    obj spectralFitting.SpectralModel
    omega {mustBeNumeric,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
end
[varargout{1:max(1,nargout)}] =...
    parAliasedSpectralDensity(obj,omega,Delta,obj.parameter);

end
function varargout = plotAcv(obj,N,Delta,varargin)
% plotAcv  Plot the autocovariance.
%   plotAcv(obj,N,Delta) will produce a plot of the autocovariance at lags 
%   determined by N and Delta.
%   Additional arguments can be provided to specify standard options for
%   plotting.
%
%   See also plot, plotSdf, plotAliasedSdf, plotEI.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}=1024
    Delta(1,1) {mustBeNumeric,mustBeNonempty} = 1
end
arguments (Repeating)
    varargin
end
if obj.useDifferencing
    M = N-1;
else
    M = N;
end
tau = (0:Delta:(M-1)*Delta)';
[varargout{1:nargout}] =...
    plot(tau,(obj.computeAutoCovariance(M,Delta)),varargin{:});
end
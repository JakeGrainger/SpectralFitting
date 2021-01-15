function varargout = plotEI(obj,N,Delta,varargin)
% plotEI  Plot the expected periodogram.
%   plotEI(obj,N,Delta) will produce a plot of the expected periodogram
%   on the decibel scale at Fourier frequencies determined by N and
%   Delta.
%   Additional arguments can be provided to specify standard options for
%   plotting.
%
%   See also plot, plotdBEI, plotSdf, plotAliasedSdf, plotAcv.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty} = 1024
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
omega = 2*pi/M/Delta*(0:floor(M/2)+1)';
[varargout{1:nargout}] = plot(omega,(obj.computeExpectedPeriodogram(M,Delta)),varargin{:});
end
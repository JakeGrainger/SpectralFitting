function varargout = plotCovI(obj,N,Delta,varargin)
% plotCovI  Plot covariance of the periodogram at Fourier frequencies.
%   plotCovI(obj,N,Delta) will produce a plot of the correlation of the
%   periodogram at Fourier frequencies determined by N and Delta.
%   Additional arguments can be provided to specify standard options for
%   plotting.
%
%   See also plot, plotCorrI, computeCovPeriodogram.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
end
arguments (Repeating)
    varargin
end
if obj.useDifferencing
    M = N-1;
else
    M = N;
end
covI = parCovPeriodogram(obj,M,Delta,obj.parameter);
covIused = covI(1:floor(M/2)+1,1:floor(M/2)+1);
omega = 2*pi/Delta/M*(0:floor(M/2))';
[varargout{1:nargout}] = imagesc(omega,omega,...
    covIused,varargin{:});
end
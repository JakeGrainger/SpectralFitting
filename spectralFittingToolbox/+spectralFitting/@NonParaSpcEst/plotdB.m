function varargout = plotdB(obj,varargin)
% method to plot objects of class NonParaSpcEst on decibel scale
[varargout{1:nargout}] = plot(obj.omega,10*log10(obj.I),varargin{:});
end
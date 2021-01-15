function varargout = plot(obj,varargin)
% overloaded plotting function of objects of class NonParaSpcEst
[varargout{1:nargout}] = plot(obj.omega,obj.I,varargin{:});
end
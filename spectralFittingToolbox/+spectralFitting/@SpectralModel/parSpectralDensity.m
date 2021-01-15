function varargout = parSpectralDensity(obj,omega,Delta,parameter)
% Computes the spectral density at omega, will compute derivatives if
% requested, for any parameter value
if obj.useDifferencing
    [Density{1:nargout}] = obj.spectralDensity(parameter,omega);
    varargout = cellfun(@(x) x.*(4*sin(omega*Delta/2).^2),Density,...
        "UniformOutput",false); % difference if necessary
else
    [varargout{1:nargout}] = obj.spectralDensity(parameter,omega);
end

end
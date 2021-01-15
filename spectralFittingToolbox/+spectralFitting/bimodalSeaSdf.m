function varargout = bimodalSeaSdf(parameter,omega)
% bimodalSeaSdf  function to evaluate the sdf for a mackay type spectra
%    [sdf,d_sdf,d2_sdf] = bimodalSeaSdf(parameter,omega) uses the
%    additiveSdf function to compute the sdf and derivatives of a process
%    defined by adding two processes with generalised JONSWAP sdfs
%    together. The vector parameter specifies the parameter of both
%    processes. The first four specify the parameters of the first
%    generalised JONSWAP in the same order as genJONSWAPsdf. The last four
%    specify the other.
%
%    See also genJONSWAPsdf, additiveSdf.
[varargout{1:nargout}] = spectralFitting.additiveSdf(parameter(1:4),@spectralFitting.genJONSWAPsdf,...
    parameter(5:end),@spectralFitting.genJONSWAPsdf,omega);
end
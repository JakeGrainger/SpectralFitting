function [sdf,d_sdf,d2_sdf] = additiveSdf(parameter1,sdf1,parameter2,sdf2,omega)
% additiveSdf evaluates the spectral density of the addition of two process
%   [sdf,d_sdf,d2_sdf] = additiveSdf(parameter1,sdf1,parameter2,sdf2,omega) 
%   evaluates the spectral density of the process X=X_1+X_2 where X_1 and
%   X_2 are two independed processes with sdf specified by function handles
%   sdf1 and sdf2 respectively. 
%   Note that the first derivative are in the same format as genJONSWAPsdf.
%   However, the second derivatives are only provided for the lower
%   triangle when they are not guarenteed to be 0. i.e. if the second
%   derivative is wrt parameters from different underlying processes, it is
%   not returned. The user must the specify this in the hessIndex property
%   of SpectralModel. See the advanced estimation tutorial for more
%   details.
%
%   See also genJONSWAPsdf, bimodalSeaSdf.
[density1{1:nargout}] = sdf1(parameter1,omega);
[density2{1:nargout}] = sdf2(parameter2,omega);
sdf = density1{1}+density2{1};
if nargout > 1
    d_sdf = [density1{2},density2{2}];
end
if nargout >2
    d2_sdf = [density1{3},density2{3}];
end
end
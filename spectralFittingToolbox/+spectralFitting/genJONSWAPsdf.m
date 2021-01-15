function [sdf,d_sdf,d2_sdf] = genJONSWAPsdf(parameter,omega)
% genJONSWAPsdf  function to evaluate a generalise JONSWAP spectra
%    [sdf,d_sdf,d2_sdf] = genJONSWAPsdf(parameter,omega) computes the
%    generalised JONSWAP spectral density with parameters defined by
%    parameter = (alpha,omega_p,gamma,r) at angular frequencies given by
%    omega. Also computes first and second derivatives, with first
%    deriviatves stacked in columns (i.e. d_sdf(:,p) is the pth
%    derivative). The second deriviative is similar, except that only the
%    lower triangle is give (as the second deriviative matrix is symetric.
%    The order given goes down eash column of the lower triangle in turn.
%
%    See also bimodalSeaSdf.
omega = abs(omega); % due to reflection
alpha = parameter(1);
omega_p = parameter(2);
gamma = parameter(3);
r = parameter(4);
s = 4;
% evaluate spectra at omega
sigma = 0.07+0.02*(omega>omega_p);
sigmasq = sigma.^2;
delta = exp(-1./(2.*sigmasq).*(omega/omega_p-1).^2);
sdf = alpha*omega.^(-r).*exp(-(r/s)*(omega/omega_p).^(-s)).*gamma.^delta;
sdf(isnan(sdf))=0; % deal with divide infinity type issues
sdf = sdf/2; % return the two sided spectral density

if nargout >= 2
    d_alpha = omega.^(-r).*exp(-(r/s)*(omega/omega_p).^(-s)).*gamma.^delta/2; d_alpha(isnan(d_alpha))=0;
    saved_omega_p = (delta.*log(gamma)./sigmasq.*omega.*(omega-omega_p)./(omega_p.^3)-omega.^(-s).*r.*omega_p.^(s-1));
    d_omega_p = sdf.*saved_omega_p;
    d_gamma = sdf.*delta/gamma;
    saved_r = (-log(omega)-(omega/omega_p).^(-s)/s);
    d_r = sdf.*saved_r;
    
    d_sdf = [d_alpha,d_omega_p,d_gamma,d_r];
    d_sdf(isnan(d_sdf))=0;
end

if nargout >= 3
    dada = zeros(size(d_alpha));
    dado = d_omega_p/alpha;
    dadg = d_gamma/alpha;
    dadr = d_r/alpha;
    dodo = d_omega_p.*saved_omega_p+sdf.*(omega.*delta*log(gamma)./(sigmasq).*(-3*omega./(omega_p.^4)+2./(omega_p.^3)+(omega-omega_p).^2.*omega./(omega_p.^6)./sigmasq)-omega_p^(s-2).*r.*(s-1)./(omega.^s));
    dodg = d_gamma.*saved_omega_p+sdf.*omega.*delta.*(omega-omega_p)/gamma/(omega_p^3)./sigmasq;
    dodr = d_r.*saved_omega_p-sdf.*(omega_p.^(s-1))./(omega.^s);
    dgdg = delta.*(d_gamma/gamma-sdf/gamma^2);
    dgdr = d_r.*delta/gamma;
    drdr = d_r.*saved_r;
    d2_sdf = [dada,dado,dadg,dadr,dodo,dodg,dodr,dgdg,dgdr,drdr];
    d2_sdf(isnan(d2_sdf))=0;
end
end
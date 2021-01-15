function varargout = simulateGaussianProcess(obj,N,Delta)
% simulateGaussianProcess  Simulate a Gaussian Process.
%   simulateGaussianProcess(obj,N,Delta) simulates a Gaussian process with
%   spectral density function defined in obj, sampled every Delta seconds,
%   recording N observations. This method uses circulant embedding.
%
%   See also randn.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
end
autocovariance = parAutoCovariance(obj,N,Delta,obj.parameter);
isPositiveSemiDefinite = false;
padding = 0; % check this
while ~isPositiveSemiDefinite
    % make s
    s = [autocovariance;zeros(padding,1);autocovariance(end-1:-1:2)]; % end-1 corresponds to lag N-2, 2 corresponds to lag 1
    % fft s
    s_tilde = fft(s);
    isPositiveSemiDefinite = min(s_tilde)>=0;
    if ~isPositiveSemiDefinite
        padding = padding + 2;
    end
    if padding > 1024
        error('failed to find an embedding');
    end
end
if padding > 0
    warning('minimal embedding not valid, %d extra values used in s',padding)
end
% gen epsilon
epsilon = randn(length(s),1)+1i*randn(length(s),1);
% make e_tilde = epsilon*sqrt(fft s/2M)
e_tilde = epsilon.*sqrt(s_tilde/length(s));
% make e = fft(e_tilde)
e = fft(e_tilde);
% subsample record
varargout{1} = real(e(1:N));
varargout{2} = imag(e(1:N));
end

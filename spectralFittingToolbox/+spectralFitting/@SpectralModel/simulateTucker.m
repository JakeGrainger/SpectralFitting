function series = simulateTucker(obj,N,Delta,M,useRandAmplitudes)
% simulateTucker  Simulate using Tuckers method.
%   simulateTucker(obj,N,Delta,M,useRandAmplitudes) will simulate an
%   approximate Gaussian process using Tuckers method. N and Delta are the
%   desired number of observations and sampling interval respectively.
%   M is the number of bins used in the frequency domain. It must be
%   bigger than N. useRandAmplitudes allows the user to choose if random
%   amplitudes are used.
%
%   See also simulateGaussianProcess.
arguments
    obj spectralFitting.SpectralModel
    N(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty}
    Delta(1,1) {mustBeNumeric,mustBeNonempty}
    M(1,1) {mustBeNumeric,mustBePositive,mustBeInteger,mustBeNonempty} = N
    useRandAmplitudes(1,1) logical = true
end
if M < N
    error("M must be larger than N.")
end
omega = 2*pi/Delta*(0:M/2)/M;
spectral_density = obj.computeSpectralDensity(omega,Delta);
sqrt_spectral_density = sqrt(spectral_density);
if useRandAmplitudes
    W = randn(1,length(spectral_density)-2);
    V = randn(1,length(spectral_density)-2);
    Zomega = (W+sqrt(-1)*V)/sqrt(2).*sqrt_spectral_density(2:end-1);
    Znyq =  randn(1,1)*sqrt_spectral_density(end);
    Z0 = randn(1,1)*sqrt_spectral_density(1);
    Z = [Z0,Zomega,Znyq,conj(Zomega(end:-1:1))];
else
    Phi = rand(1,length(spectral_density))*2*pi;
    Zplus = sqrt_spectral_density.*exp(1i*Phi);
    Z = [Zplus,conj(Zplus(end-1:-1:2))];
end
overSampledTimeseries = real(sqrt(2*pi/M/Delta)*fft(Z));
series = overSampledTimeseries(1:N);

end

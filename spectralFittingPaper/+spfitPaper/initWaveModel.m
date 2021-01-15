function initialParameter = initWaveModel(record,Delta)
% initWaveModel  Choose initial values for parameters.
%   initWaveModel(record,Delta) produces and initial guess for the
%   parameters of a generalised JONSWAP based on the provided record.
%   Useful helper function for the experiments in paper.
%% get the peak
Iobs = spectralFitting.Periodogram(record,Delta);
[~,maxIndex] = max(Iobs.I);
omega_p = Iobs.omega(maxIndex);
%% get r
r = -log(Iobs.omega(maxIndex+floor((end-maxIndex)/2):end)) \ log(Iobs.I(maxIndex+floor((end-maxIndex)/2):end));
%% get gamma
gamma = 3;
%% get alpha
waveModelAlphaOne = spectralFitting.SpectralModel(@spectralFitting.genJONSWAPsdf,...
    [1,omega_p,gamma,r],[0;0;1;1],Inf(4,1));
cov = waveModelAlphaOne.computeAutoCovariance(length(record),Delta);
alpha = var(record)/cov(1); % i.e. obsVar = modelVar and model variance is alpha times variance when alpha is set to 1.
%% make vector
initialParameter = [alpha,omega_p,gamma,r];
end
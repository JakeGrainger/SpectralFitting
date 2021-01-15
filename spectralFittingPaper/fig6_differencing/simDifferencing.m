%%% script to simulate the data used in to produce figure 2 and 3 from the
%%% paper.
clear
%% set up random number generator
rng(1192,'twister')              % this seed will reproduce the figure from
                                 % the paper, it was selected randomly. Feel 
                                 % free to change it and try other values
generatorSettings = rng;         % saves random number generator settings

%% set parameters
nReps = 1000;                     % number of repetitions
N = 2048*4;                      % number of observations per record
Delta = 1/4;                     % sampling interval
trueParameter = [0.7;0.8;2.5;5]; % true parameter vector

%% set up model
waveSpecModel = spectralFitting.SpectralModel(...
    @spectralFitting.genJONSWAPsdf,trueParameter,[0;0;1;1],[10;pi/Delta;20;10]);

%% simulate series and save
series4Hz = nan(N,nReps);
for i = 2:2:nReps
    [series4Hz(:,i-1),series4Hz(:,i)] = waveSpecModel.simulateGaussianProcess(N,Delta);
end
series1Hz = series4Hz(1:4:end,:);
save('datDifferencing.mat')
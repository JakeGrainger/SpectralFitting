%%% script to simulate the data used in to produce figure 4 from the paper.
clear
%% set up random number generator
rng(5547,'twister')              % this seed will reproduce the figure from
                                 % the paper, it was selected randomly. Feel 
                                 % free to change it and try other values 
                                 % (this script will take a while to run as 
                                 % the maximum likelihood is slow, so use 
                                 % with caution - ideally on a cluster.
                                 % (Change for to parfor).
generatorSettings = rng;         % saves random number generator settings

%% set parameters
nReps = 1000;                     % number of repetitions
N = 13824;                       % number of observations per record
Delta = 1/1.28;                  % sampling interval
trueParameter = [0.7;0.7;3.3;4]; % true parameter vector

%% set up model
waveSpecModel = spectralFitting.SpectralModel(...
    @spectralFitting.genJONSWAPsdf,trueParameter,[0;0;1;1],[10;pi/Delta;20;10]);

%% simulate series and save
seriesStore = nan(N,nReps);
for iRep = 2:2:nReps
    [seriesStore(:,iRep-1),seriesStore(:,iRep)] =...
        waveSpecModel.simulateGaussianProcess(N,Delta);
end
save('datThreeHour.mat')
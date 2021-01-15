%% Hyperparameters
N = 2304;
Delta = 1/1.28;
nReps = 1000;
%% setup parameter choices
alphaVals = 0.7; % the choice of alpha wont make a difference as it just does scaling.
omegapVals = [0.7,0.9,1.2];
gammaVals = [1,2,3.3,5];
rVals = [4,5];
%% build cell of parameter vectors
count = 1;
nTruePar = length(omegapVals)*length(gammaVals)*length(rVals);
parameterChoice = cell(nTruePar,1);
for ii = 1:length(omegapVals)
    for jj = 1:length(gammaVals)
        for kk = 1:length(rVals)
            parameterChoice{count} = [alphaVals;omegapVals(ii);gammaVals(jj);rVals(kk)];
            count = count + 1;
        end
    end
end
%%
waveSpecModel = spectralFitting.SpectralModel(...
    @spectralFitting.genJONSWAPsdf,parameterChoice{1},...
    [0;0;1;1],[10;2;30;15]);

%%
figure();clf;
for ii = 1:nTruePar
    waveSpecModel.parameter = parameterChoice{ii};
    subplot(6,4,ii); waveSpecModel.plotCorrI(N, Delta);
    title(sprintf("a = %0.1f, om = %0.1f, ga = %0.1f, r = %0.1f", parameterChoice{ii}))
end
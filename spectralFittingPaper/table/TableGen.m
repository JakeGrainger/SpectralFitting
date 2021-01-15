clear
rng(5367)
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
% set low freq cutoff
count = 1;
lowFreqCutoff = cell(nTruePar,1);
for ii = 1:length(omegapVals)
    for jj = 1:length(gammaVals)
        for kk = 1:length(rVals)
            if parameterChoice{count}(2) == 0.7
                lowFreqCutoff{count} = 0.4;
            elseif parameterChoice{count}(2) == 0.9
                lowFreqCutoff{count} = 0.6;
            elseif parameterChoice{count}(2) == 1.2
                lowFreqCutoff{count} = 0.7;
            else
                error('Failed to allocate low freq cutoff')
            end
            count = count + 1;
        end
    end
end
%%
% loop over parameter choice
fprintf("Running %d iterations of fitting procedure.\n",nTruePar)
fullFits = nan(4, nReps, 3, nTruePar);
for iTruePar = 1:nTruePar
    tic
    fullFits(:,:,:,iTruePar) = simFitPlotExperiment(parameterChoice{iTruePar},N,Delta,nReps,iTruePar, lowFreqCutoff{iTruePar});
    itrTime = toc;
    fprintf("Finished iteration %d in %0.2f seconds.\n",iTruePar,itrTime)
end

%% save
save('LargeExperiment.mat')

%% function to perform simulation, fitting and plotting
function fittedValues = simFitPlotExperiment(trueParameter,N,Delta,nReps,expNum, lfCutoff)
% set up model
waveSpecModel = spectralFitting.SpectralModel(...
    @spectralFitting.genJONSWAPsdf,trueParameter,...
    [0;0.4;1;1],[10;2;20;10]);
waveSpecModel.bartlettWindowSize = 100/Delta; % set Bartlett window size
fitIndex = 2*pi/Delta/N*(1:N/2+1)' > lfCutoff; % remove low freq from fit
fitIndexBart = 2*pi/Delta/waveSpecModel.bartlettWindowSize*...
    (1:waveSpecModel.bartlettWindowSize/2+1)' > lfCutoff; % remove for bartlett
waveSpecModel.fitRoutine = 'fmincon'; % faster than fminsearch for applicable methods (when we have derivatives and hessian)
% simulate series and save
seriesStore = nan(N,nReps);
for iRep = 2:2:nReps
    [seriesStore(:,iRep-1),seriesStore(:,iRep)] =...
        waveSpecModel.simulateGaussianProcess(N,Delta);
end
% perform fitting
fittingMethods = {'leastSquares', 'bartlettLeastSquares',...
    'debiasedWhittle'};
fittedValues = nan(length(trueParameter),nReps,length(fittingMethods));
for iFit = 1:length(fittingMethods)
    switch fittingMethods{iFit}
        case 'bartlettLeastSquares'
            fitIndexUsed = fitIndexBart;
        otherwise
            fitIndexUsed = fitIndex;
    end
    waveSpecModel.objectiveMethod = fittingMethods{iFit};
    for jRep = 1:nReps
        % choose initial guess
        waveSpecModel.parameter = altInitWaveModel(seriesStore(:,jRep),Delta);
        % fit using initial guess
        fittedValues(:,jRep,iFit) = waveSpecModel.estimateParameter(seriesStore(:,jRep),Delta,fitIndexUsed);
        if mod(jRep, 100) == 0
            fprintf("Done iteration %d of fitting %d of experiment %d\n", jRep, iFit, expNum)
        end
    end
end
end

% in practice, fixing an initial r = 4.5 can be better
function initialParameter = altInitWaveModel(record,Delta)
% initWaveModel  Choose initial values for parameters.
%   initWaveModel(record,Delta) produces and initial guess for the
%   parameters of a generalised JONSWAP based on the provided record.
%   Useful helper function for the experiments in paper.
%% get the peak
Iobs = spectralFitting.Periodogram(record,Delta);
[~,maxIndex] = max(Iobs.I);
omega_p = Iobs.omega(maxIndex);
%% get r
r = 4.5;
%% get gamma
gamma = 3;
%% get alpha
waveModelAlphaOne = spectralFitting.SpectralModel(@spectralFitting.genJONSWAPsdf,...
    [1,omega_p,gamma,r],[0;0;1;1],Inf(4,1));
cov = waveModelAlphaOne.computeAutoCovariance(length(record),Delta);
alpha = var(record)/cov(1); % i.e. obsVar = modelVar and model variance is alpha times variance when alpha is set to 1.
%% make vector
initialParameter = [alpha;omega_p;gamma;r];
end
%%% script to fit a generalised JONSWAP to produce figure 5 from the paper.
clear
totalStart = tic;
%% load simulated series
load('datGammaEqualOne')

%% set fitting settings
waveSpecModel.bartlettWindowSize = 100/Delta; % set Bartlett window size
fitIndex = 2*pi/Delta/N*(1:N/2+1)' > 0.4; % remove low freq from fit
fitIndexBart = 2*pi/Delta/waveSpecModel.bartlettWindowSize*...
    (1:waveSpecModel.bartlettWindowSize/2+1)' > 0.4; % remove for bartlett

%% perform fitting
fittingMethods = {'leastSquares', 'bartlettLeastSquares',...
    'debiasedWhittle'};
fittingSpeed = nan(nReps,length(fittingMethods));
fittedValues = nan(length(trueParameter),nReps,length(fittingMethods));
for iFit = 1:length(fittingMethods)
    switch fittingMethods{iFit}
        case 'bartlettLeastSquares'
            fitIndexUsed = fitIndexBart;
        otherwise
            fitIndexUsed = fitIndex;
    end
    waveSpecModel.objectiveMethod = fittingMethods{iFit};
    waveSpecModel.fitRoutine = 'fmincon';
    for jRep = 1:nReps
        % choose initial guess
        waveSpecModel.parameter = altInitWaveModel(seriesStore(:,jRep),Delta);
        % fit using initial guess
        thisTime = tic;
        fittedValues(:,jRep,iFit) = waveSpecModel.estimateParameter(seriesStore(:,jRep),Delta,fitIndexUsed);
        % save time taken
        fittingSpeed(jRep,iFit) = toc(thisTime);
        % display iteration
        disp(jRep)
    end
    
end

%% save
save('fittedGammaEqualOne.mat')
toc(totalStart)

%% alternate initialisation function that we found often made optimisation faster as we know in practice r is close to 4.5, this is often better than the log version as that tends to over estimate r (though it doesnt really make a huge difference what staring values are used)
function initialParameter = altInitWaveModel(record,Delta)

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
initialParameter = [alpha,omega_p,gamma,r];
end
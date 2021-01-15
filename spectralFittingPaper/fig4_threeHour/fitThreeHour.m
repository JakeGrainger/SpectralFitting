%%% script to fit a generalised JONSWAP to produce figure 4 from the paper.
clear
totalStart = tic;
%% load simulated series
load('datThreeHour')

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
        parWaveSpecModel = waveSpecModel;
        parWaveSpecModel.parameter = spfitPaper.initWaveModel(seriesStore(:,jRep),Delta);
        % fit using initial guess
        thisTime = tic;
        fittedValues(:,jRep,iFit) = parWaveSpecModel.estimateParameter(seriesStore(:,jRep),Delta,fitIndexUsed);
        % save time taken
        fittingSpeed(jRep,iFit) = toc(thisTime);
        % print iteration number
        disp(jRep)
    end
    
end

%% save
clear seriesStore % series are too large to store on GitHub
save('fittedThreeHour.mat')
toc(totalStart)
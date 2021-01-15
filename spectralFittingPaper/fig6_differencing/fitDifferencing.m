%%% script to fit a generalised JONSWAP to produce figure 2 and 3 from the
%%% paper.
clear
tic
%% load simulated series
load('datDifferencing')

%% fit model
sampleInt = {1,1/4};
series = {series1Hz,series4Hz};
dwValues = cell(length(sampleInt),1); lsValues = cell(length(sampleInt),1); blsValues = cell(length(sampleInt),1);
%% normal
waveSpecModel.fitRoutine = 'fmincon';
waveSpecModel.useDifferencing = false;
for ii = 1:length(sampleInt)
    M = size(series{ii},1);
    fitIndex = 2*pi/sampleInt{ii}/M*(1:M/2+1)' > 0.4;
    waveSpecModel.objectiveMethod = 'debiasedWhittle';
    dwValues{ii} = waveSpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    waveSpecModel.objectiveMethod = 'leastSquares';
    lsValues{ii} = waveSpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    fitIndex = 2*pi/sampleInt{ii}/(100/sampleInt{ii})*(1:(100/sampleInt{ii})/2+1)' > 0.4;
    waveSpecModel.bartlettWindowSize = 100/sampleInt{ii};
    waveSpecModel.objectiveMethod = 'bartlettLeastSquares';
    blsValues{ii} = waveSpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
end

%% differenced
dwDiff = cell(length(sampleInt),1); lsDiff = cell(length(sampleInt),1); blsDiff = cell(length(sampleInt),1);
waveSpecModel.useDifferencing = true;
for ii = 1:length(sampleInt)
    M = size(series{ii},1)-1;
    fitIndex = 2*pi/sampleInt{ii}/M*(1:M/2+1)' > 0.5;
    waveSpecModel.objectiveMethod = 'debiasedWhittle';
    dwDiff{ii} = waveSpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    waveSpecModel.objectiveMethod = 'leastSquares';
    lsDiff{ii} = waveSpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    waveSpecModel.objectiveMethod = 'bartlettLeastSquares';
    fitIndex = 2*pi/sampleInt{ii}/(100/sampleInt{ii})*(1:(100/sampleInt{ii})/2+1)' > 0.5;
    waveSpecModel.bartlettWindowSize = 100/sampleInt{ii};
    blsDiff{ii} = waveSpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
end
toc

%% save
clear series series1Hz series4Hz
save('fittedDifferencing.mat')

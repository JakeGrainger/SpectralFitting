clear
load('LargeExperiment.mat')
%% calculate bias, sd and rmse
biasSep = nan(4, 3, nTruePar);
stdSep = nan(4, 3, nTruePar);
rmseSep = nan(4, 3, nTruePar);
for iTruePar = 1:nTruePar
    fittedValues = squeeze(fullFits(:, :, :, iTruePar));
    trueParameter = parameterChoice{iTruePar};
    biasSep(:, :, iTruePar) = squeeze(abs(mean(fittedValues - trueParameter, 2)))./trueParameter;
    stdSep(:, :, iTruePar) = squeeze(std(fittedValues - trueParameter, [], 2))./trueParameter;
    rmseSep(:, :, iTruePar) = squeeze(sqrt(mean((fittedValues - trueParameter).^2, 2)))./trueParameter;
end
biasTot = mean(biasSep, 3);
stdTot = mean(stdSep, 3);
rmseTot = mean(rmseSep, 3);

%%
format bank
disp('bias')
disp(biasTot*100)

%%
disp('std')
disp(stdTot*100)

%%
disp('rmse')
disp(rmseTot*100)

%% average over pars
disp('bias average')
disp(mean(biasTot*100,1))

%%
disp('std average')
disp(mean(stdTot*100,1))

%%
disp('rmse average')
disp(mean(rmseTot*100,1))

format short

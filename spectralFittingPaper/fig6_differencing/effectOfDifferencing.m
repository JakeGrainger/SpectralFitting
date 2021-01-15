%% get differenced and non differenced estimates for both LS and DW
rng(3190,'twister')
% settings
tic
nReps = 500;
N = 4096;
Delta = 1/4;
trueParameter = [0.7;0.7;3.3;5];
% make model
mackaySpecModel = spectralFitting.SpectralModel(@spectralFitting.genJONSWAPsdf,...
    trueParameter,[0;0;1;1],[10;pi/Delta;30;15]);
mackaySpecModel.fitRoutine = 'fmincon';
% do simulations
series4Hz = nan(N,nReps);
for i = 2:2:nReps
    [series4Hz(:,i-1),series4Hz(:,i)] = mackaySpecModel.simulateGaussianProcess(N,Delta);
end
series1Hz = series4Hz(1:4:end,:);
% series2Hz = series4Hz(1:2:end,:);
disp("Finished simulation.")
%% fit model
sampleInt = {1,1/4};
series = {series1Hz,series4Hz};
dwValues = cell(length(sampleInt),1); lsValues = cell(length(sampleInt),1); blsValues = cell(length(sampleInt),1);
% normal
mackaySpecModel.useDifferencing = false;
for ii = 1:length(sampleInt)
    M = size(series{ii},1);
    fitIndex = 2*pi/sampleInt{ii}/M*(1:M/2+1)' > 0.38;
    mackaySpecModel.objectiveMethod = 'debiasedWhittle';
    dwValues{ii} = mackaySpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    mackaySpecModel.objectiveMethod = 'leastSquares';
    lsValues{ii} = mackaySpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    fitIndex = 2*pi/sampleInt{ii}/(128/sampleInt{ii})*(1:(128/sampleInt{ii})/2+1)' > 0.38;
    mackaySpecModel.bartlettWindowSize = 128/sampleInt{ii};
    mackaySpecModel.objectiveMethod = 'bartlettLeastSquares';
    blsValues{ii} = mackaySpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
end
disp("finished normal fits.")
% differenced
dwDiff = cell(length(sampleInt),1); lsDiff = cell(length(sampleInt),1); blsDiff = cell(length(sampleInt),1);
mackaySpecModel.useDifferencing = true;
for ii = 1:length(sampleInt)
    M = size(series{ii},1)-1;
    fitIndex = 2*pi/sampleInt{ii}/M*(1:M/2+1)' > 0.38;
    mackaySpecModel.objectiveMethod = 'debiasedWhittle';
    dwDiff{ii} = mackaySpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    mackaySpecModel.objectiveMethod = 'leastSquares';
    lsDiff{ii} = mackaySpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
    mackaySpecModel.objectiveMethod = 'bartlettLeastSquares';
    fitIndex = 2*pi/sampleInt{ii}/(128/sampleInt{ii})*(1:(128/sampleInt{ii})/2+1)' > 0.38;
    mackaySpecModel.bartlettWindowSize = 128/sampleInt{ii};
    blsDiff{ii} = mackaySpecModel.estimateParameter(series{ii},sampleInt{ii},fitIndex);
end
disp("finished differenced fits.")
toc
%% fig settings
opts = spfitPaper.appliedPlotSettings();
%%
for ii = 1:2
clear fittedValues
fittedValues(:,:,1) = lsValues{ii};
fittedValues(:,:,2) = blsValues{ii};
fittedValues(:,:,3) = dwValues{ii};
fittedValues(:,:,4) = lsDiff{ii};
fittedValues(:,:,5) = blsDiff{ii};
fittedValues(:,:,6) = dwDiff{ii};
%
parameterName = {'\alpha','\omega_p','\gamma','r'};
plotOptions = spfitPaper.appliedBoxSettings([1,2,5]);
plotOptions{6} = repmat(plotOptions{6},2,1);
plotOptions{4} = 1:6;
fig = figure();clf;hold on;
ax = spfitPaper.multipleBox(permute(fittedValues,[2,3,1]),{'LS','BLS','DW','LSd','BLSd','DWd'},plotOptions{:});
for jPar = 1:4
    ax(jPar).YLabel.String = parameterName{jPar};
    yline(ax(jPar),trueParameter(jPar),'--');
    grid(ax(jPar),'on');
    if jPar == 1
        set(ax(jPar),'ylim',[0.1,1.3]);
    elseif jPar == 2
        set(ax(jPar),'ylim',[0.63,0.79]);
    elseif jPar > 2
        set(ax(jPar),'ylim',[0.5,10.5]);
    end
end
fig.Units               = 'centimeters';
fig.Position(3)         = opts.size/2;
fig.Position(4)         = 16;
set(fig.Children, ...
    'FontName',     opts.font, ...
    'FontSize',     opts.font_size);
set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
fig.PaperPositionMode   = 'auto';
print(sprintf('differencingBox%d',1/sampleInt{ii}), '-dpng', '-r600' )
end
% close all

%%
BartStore = nan(128*2+1,nReps);
PerStore = nan(2049,nReps);
BartDiff = nan(128*2+1,nReps);
PerDiff = nan(2048,nReps);
for ii = 1:nReps
    BartStore(:,ii) = spectralFitting.BartlettPeriodogram(series4Hz(:,ii),1/4,128*4).I;
    PerStore(:,ii) = spectralFitting.Periodogram(series4Hz(:,ii),1/4).I;
    BartDiff(:,ii) = spectralFitting.BartlettPeriodogram(diff(series4Hz(:,ii)),1/4,128*4).I;
    PerDiff(:,ii) = spectralFitting.Periodogram(diff(series4Hz(:,ii)),1/4).I;
end
Delta = 1/4;
genStore = {BartStore,PerStore,BartDiff,PerDiff};
plotName = {'Bn','Pn','Bd','Pd'};
for ii = 1:4
        fig = figure();clf;
        colormap(cool)
        imagesc([0,pi/Delta],[0,pi/Delta],corrcov(cov(genStore{ii}')))
        xlim = ([0,pi/Delta]); ylim = ([0,pi/Delta]);
        xticks((0:2:4)*pi)
        xticklabels(["0","2\pi","4\pi"])
        yticks((0:2:4)*pi)
        yticklabels(["0"," 2\pi"," 4\pi"])
        grid on; box on
        fig.Units               = 'centimeters';
        fig.Position(3)         = 7.75;
        fig.Position(4)         = 4.5;
        if ii == 4
            colorbar('location','westoutside','visible','off')
            colorbar('location','eastoutside','Limits',[0,1])
        else
            colorbar('location','westoutside','visible','off')
            colorbar('location','eastoutside','visible','off')
        end
        set(fig.Children, ...
            'FontName',     opts.font, ...
            'FontSize',     opts.font_size);
        set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
        fname = sprintf('corr_%s',plotName{ii});
        print(fname, '-dpng', '-r600' )
end

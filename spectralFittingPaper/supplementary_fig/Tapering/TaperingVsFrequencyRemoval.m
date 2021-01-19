clear
%% set parameters
nReps = 1000;                     % number of repetitions
N = 2304;                        % number of observations per record
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
%%
filts=dpss(N,4)'; hh1=filts(1,:)';
C = abs(fft(seriesStore ,[], 1)).^2/N/2/pi*Delta; C = C';
D = abs(fft(seriesStore.*hh1 ,[], 1)).^2/2/pi*Delta; D = D';
filts2=dpss(N-1,4)'; hh2=filts2(1,:)';
E = abs(fft(diff(seriesStore,1,1) ,[], 1)).^2/N/2/pi*Delta; E = E';
F = abs(fft(diff(seriesStore,1,1).*hh2 ,[], 1)).^2/2/pi*Delta; F = F';
P = spectralFitting.Periodogram(seriesStore(:,1), Delta);
fig = figure();clf;
subplot(1,4,1); imagesc(P.omega,P.omega,corr(C(:, 1:end/2))); title("Periodogram")
subplot(1,4,2); imagesc(P.omega,P.omega,corr(D(:, 1:end/2))); title("Tapered Periodogram")
subplot(1,4,3); imagesc(P.omega,P.omega,corr(E(:, 1:floor(end/2)))); title("Periodogram Differenced")
subplot(1,4,4); imagesc(P.omega,P.omega,corr(F(:, 1:floor(end/2)))); title("Tapered Periodogram Differenced")
fig.Units               = 'centimeters';
fig.Position(3)         = 25;
fig.Position(4)         = 5;
print -dpng -r600 tapering_correlation_periodogram
%%
fitIndex = 2*pi/Delta/N*(1:N/2+1)' > 0.4;
parEst = nan(4, size(seriesStore, 2), 3);
waveSpecModel.useDifferencing = false;
tic
waveSpecModel.fitRoutine = 'fmincon';
waveSpecModel.objectiveMethod = 'debiasedWhittle';
parEst(:,:,1) = waveSpecModel.estimateParameter(seriesStore, Delta);
toc

tic
parEst(:,:,2) = waveSpecModel.estimateParameter(seriesStore, Delta, fitIndex);
toc

tic
waveSpecModel.objectiveMethod = 'taperedDW';
parEst(:,:,3) = waveSpecModel.estimateParameter(seriesStore, Delta);
toc

%%
M = N-1;
diffFitIndex = 2*pi/Delta/M*(1:floor(M/2)+1)' > 0.4;
diffParEst = nan(4, size(seriesStore, 2), 3);
waveSpecModel.useDifferencing = true;

tic
waveSpecModel.fitRoutine = 'fmincon';
waveSpecModel.objectiveMethod = 'debiasedWhittle';
diffParEst(:,:,1) = waveSpecModel.estimateParameter(seriesStore, Delta);
toc

tic
diffParEst(:,:,2) = waveSpecModel.estimateParameter(seriesStore, Delta, diffFitIndex);
toc

tic
waveSpecModel.objectiveMethod = 'taperedDW';
diffParEst(:,:,3) = waveSpecModel.estimateParameter(seriesStore, Delta);
toc



%%
% for iPar = 1:4
%     subplot(4,1,iPar); boxplot(squeeze(parEst(iPar,:,:)), ["DW", "DWremoved", "DWtapered"])
%     ylabel(parName(iPar));
%     if iPar == 3
%         ylim([0,10])
%     end
%     yline(trueParameter(iPar))
% end
%% get options for plotting
boxOptions = spfitPaper.appliedBoxSettings(1:6);
opts = spfitPaper.appliedPlotSettings();

%% plot figure
fullParEst = nan(4, size(seriesStore, 2), 6);
fullParEst(:, :, 1:3) = parEst;
fullParEst(:, :, 4:6) = diffParEst;

fittingNames = ["DW", "DWremoved", "DWtapered", "differenced DW", "differenced DWremoved", "differenced DWtapered"];
parameterName = {'\alpha','\omega_p','\gamma','r'};
fig = figure();clf;hold on;
ax = spfitPaper.multipleBox(permute(fullParEst(:,:,:),[2,3,1]),...
    fittingNames,boxOptions{:});
for jPar = 1:4
    ax(jPar).YLabel.String = parameterName{jPar};
    yline(ax(jPar),trueParameter(jPar),'--');
    grid(ax(jPar),'on');
end
xtickangle(ax(4),30);
fig.Units               = 'centimeters';
fig.Position(3)         = opts.size;
fig.Position(4)         = 20;
set(fig.Children, ...
    'FontName',     opts.font, ...
    'FontSize',     opts.font_size);
set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
fig.PaperPositionMode   = 'auto';
print -dpng -r600 tapering_effects

clear
%%

sdf = @spectralFitting.genJONSWAPsdf;
theta = [0.7;0.7;3.3;4];
lb = [0;0.4;1;1];
Delta = 1/1.28;
ub = [10;pi/Delta;20;10];
waveSpecModel = spectralFitting.SpectralModel(sdf,theta,lb,ub);
N = 2304;
nReps = 1000;

omega = 2*pi/Delta/N*(1:floor(N/2)+1)';
fitIndex = omega > 0.4;

fitIndexBart = 2*pi/Delta/waveSpecModel.bartlettWindowSize*...
    (1:waveSpecModel.bartlettWindowSize/2+1)' > 0.4;
%%
parameterEstimate = nan(nReps, 4, 3);
waveSpecModel.fitRoutine = 'fmincon';
acv = computeAutoCovariance(waveSpecModel, N, Delta);
C = toeplitz(acv);
for ii = 1:nReps
    df = 16;
    ts = mvtrnd(C, df)* sqrt((df-2)/df*acv(1));
    waveSpecModel.objectiveMethod = 'leastSquares';
    parameterEstimate(ii,:,1) = waveSpecModel.estimateParameter(ts', Delta, fitIndex);
    waveSpecModel.objectiveMethod = 'bartlettLeastSquares';
    parameterEstimate(ii,:,2) = waveSpecModel.estimateParameter(ts', Delta, fitIndexBart);
    waveSpecModel.objectiveMethod = 'debiasedWhittle';
    parameterEstimate(ii,:,3) = waveSpecModel.estimateParameter(ts', Delta, fitIndex);
    if mod(ii, 100) == 0
        fprintf("Finished iteration %d\n", ii)
    end
end
%% get options for plotting
boxOptions = spfitPaper.appliedBoxSettings([1,2,5]);
opts = spfitPaper.appliedPlotSettings();
nPar = 4;
%% plot figure
parameterName = {'\alpha','\omega_p','\gamma','r'};
fittingName = {'Least Squares', 'Bartlett LS', 'Debiased Whittle'};
fig = figure();clf;hold on;
ax = spfitPaper.multipleBox(permute(parameterEstimate,[1,3,2]),...
    fittingName,boxOptions{:});
for jPar = 1:nPar
    ax(jPar).YLabel.String = parameterName{jPar};
    yline(ax(jPar),theta(jPar),'--');
    grid(ax(jPar),'on');
end
fig.Units               = 'centimeters';
fig.Position(3)         = opts.size;
fig.Position(4)         = 20;
set(fig.Children, ...
    'FontName',     opts.font, ...
    'FontSize',     opts.font_size);
set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
fig.PaperPositionMode   = 'auto';
print -dpng -r600 t_fits_df_16
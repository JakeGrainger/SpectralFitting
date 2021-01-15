%%% script to simulate the data used in to produce figure 5 from the paper.
clear
%% set up random number generator
rng(7281,'twister')              % this seed will reproduce the figure from
                                 % the paper, it was selected randomly. Feel 
                                 % free to change it and try other values 
                                 % (this script will take a while to run as 
                                 % the maximum likelihood is slow, so use 
                                 % with caution - ideally on a cluster.
                                 % (Change for to parfor).
generatorSettings = rng;         % saves random number generator settings

%% set parameters
nReps = 1000;                     % number of repetitions
N = 2304;                        % number of observations per record
Delta = 1/1.28;                  % sampling interval
trueParameter = [0.7;0.7;3.3;4]; % true parameter vector

%% set up model
waveSpecModel = spectralFitting.SpectralModel(...
    @spectralFitting.genJONSWAPsdf,trueParameter,[0;0;1;1],[10;pi/Delta;20;10]);
%% simulate
seriesStore = nan(N,nReps);
for iRep = 2:2:nReps
    [seriesStore(:,iRep-1),seriesStore(:,iRep)] =...
        waveSpecModel.simulateGaussianProcess(N,Delta);
end
%% fit
%% set fitting settings
waveSpecModel.bartlettWindowSize = 100/Delta; % set Bartlett window size
omega = 2*pi/Delta/N*(1:N/2+1)';
fitIndex = omega > 0.4; % remove low freq from fit
bartOmega = 2*pi/Delta/waveSpecModel.bartlettWindowSize*...
    (1:waveSpecModel.bartlettWindowSize/2+1)';
fitIndexBart = bartOmega > 0.4; % remove for bartlett

%% perform fitting
fittingMethods = {'leastSquares', 'bartlettLeastSquares',...
    'debiasedWhittle'};
fittedValues = nan(length(trueParameter),nReps,7);
fitInd = [1, 2, 7];
for iFit = 1:3
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
        fittedValues(:,jRep,fitInd(iFit)) = waveSpecModel.estimateParameter(seriesStore(:,jRep),Delta,fitIndexUsed);
        % display iteration
        if mod(jRep,200) == 0
            fprintf("Finised %d of part %d\n", jRep, iFit)
        end
    end 
end
for jRep = 1:nReps
    P = spectralFitting.Periodogram(seriesStore(:, jRep), Delta);
    B = spectralFitting.BartlettPeriodogram(seriesStore(:, jRep), Delta, waveSpecModel.bartlettWindowSize);
    iniParameter = altInitWaveModel(seriesStore(:,jRep),Delta);
    
    obj = @(theta) logLS(waveSpecModel, P.I(fitIndex), P.omega(fitIndex), theta, Delta);
    fittedValues(:,jRep,3) = fminsearch(obj, iniParameter, waveSpecModel.fminsearchOpts);
    
    obj = @(theta) logLS(waveSpecModel, B.I(fitIndexBart), B.omega(fitIndexBart), theta, Delta);
    fittedValues(:,jRep,4) = fminsearch(obj, iniParameter, waveSpecModel.fminsearchOpts);
    
    obj = @(theta) aliasedLogLS(waveSpecModel, P.I(fitIndex), P.omega(fitIndex), theta, Delta);
    fittedValues(:,jRep,5) = fminsearch(obj, iniParameter, waveSpecModel.fminsearchOpts);
    
    obj = @(theta) aliasedLogLS(waveSpecModel, B.I(fitIndexBart), B.omega(fitIndexBart), theta, Delta);
    fittedValues(:,jRep,6) = fminsearch(obj, iniParameter, waveSpecModel.fminsearchOpts);
    if mod(jRep, 200) == 0
        fprintf("Finised %d of long\n", jRep)
    end
end

%% plot
plotcolor = nan(7, 3);
plotcolor(1:2,:) = spfitPaper.appliedColor([1, 2]);
plotcolor(7,:) = spfitPaper.appliedColor(5);
plotcolor(3:6,:) = [0,0,255;
    0,255,0;
    255,0,0;
    139,69,19]/256;


plotOptions = {"PlotStyle","traditional",...
    "ColorGroup",1:7,...
    "Colors",plotcolor,...
    "BoxStyle","outline",...
    "MedianStyle","line",...
    "Symbol","+"};
opts = spfitPaper.appliedPlotSettings();
nPar = 4;
parameterName = {'\alpha','\omega_p','\gamma','r'};
fittingName = {'Least Squares', 'Bartlett LS', 'log Least Squares','log Bartlett LS', 'aliased log Least Squares','aliased log Bartlett LS', 'Debiased Whittle'};
fig = figure();clf;hold on;

ax = spfitPaper.multipleBox(permute(fittedValues(:,:,:),[2,3,1]),...
    fittingName,plotOptions{:});
for jPar = 1:nPar
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
print -dpng -r600 log_LS_box

%% plot without first two
plotcolor = nan(5, 3);
plotcolor(5,:) = spfitPaper.appliedColor(5);
plotcolor(1:4,:) = [0,0,255;
    0,255,0;
    255,0,0;
    139,69,19]/256;


plotOptions = {"PlotStyle","traditional",...
    "ColorGroup",1:5,...
    "Colors",plotcolor,...
    "BoxStyle","outline",...
    "MedianStyle","line",...
    "Symbol","+"};
opts = spfitPaper.appliedPlotSettings();
nPar = 4;
parameterName = {'\alpha','\omega_p','\gamma','r'};
fittingName = {'log Least Squares','log Bartlett LS', 'aliased log Least Squares','aliased log Bartlett LS', 'Debiased Whittle'};
fig = figure();clf;hold on;

ax = spfitPaper.multipleBox(permute(fittedValues(:,:,3:end),[2,3,1]),...
    fittingName,plotOptions{:});
for jPar = 1:nPar
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
print -dpng -r600 log_LS_box_noLS

%% helper functions
function val = logLS(model, I, omega, theta, Delta)
    model.parameter = theta;
    sdf = model.computeSpectralDensity(omega, Delta);
    val = sum((log(sdf) - log(I)).^2);
end
function val = aliasedLogLS(model, I, omega, theta, Delta)
    model.parameter = theta;
    sdf = model.computeAliasedSpectralDensity(omega, Delta);
    val = sum((log(sdf) - log(I)).^2);
end

function initialParameter = altInitWaveModel(record, Delta)

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
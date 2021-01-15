%%% script to plot figure 2 and 3 from the paper.
clear
%% load fitted data
load('fittedHalfHour')

%% add timings as extra parameter
nPar = length(trueParameter);
fullValues = nan(nPar+1,nReps,length(fittingMethods));
fullValues(1:nPar,:,:) = fittedValues;
fullValues(end,:,:) = log10(fittingSpeed);

%% get options for plotting
boxOptions = spfitPaper.appliedBoxSettings(1:6);
opts = spfitPaper.appliedPlotSettings();

%% plot figure 2
parameterName = {'\alpha','\omega_p','\gamma','r','log10 runtime'};
fittingName = {'Least Squares', 'Bartlett LS', 'Whittle', 'Aliased Whittle', 'Debiased Whittle', 'Maximum Likelihood'};
fig = figure();clf;hold on;
ax = spfitPaper.multipleBox(permute(fullValues(:,:,:),[2,3,1]),...
    fittingName,boxOptions{:});
for jPar = 1:nPar+1
    ax(jPar).YLabel.String = parameterName{jPar};
    if jPar <= nPar
        yline(ax(jPar),trueParameter(jPar),'--');
    end
    grid(ax(jPar),'on');
end
xtickangle(ax(5),30);
fig.Units               = 'centimeters';
fig.Position(3)         = opts.size;
fig.Position(4)         = 20;
set(fig.Children, ...
    'FontName',     opts.font, ...
    'FontSize',     opts.font_size);
set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
fig.PaperPositionMode   = 'auto';
print -dpng -r600 half_hour_box

%% plot figure 3
fig = figure();clf;
[S,ax,~,H] = plotmatrix(squeeze(fittedValues(:,:,5))');
for iPar = 1:nPar
    ax(iPar,1).YLabel.String=parameterName{iPar};
    ax(4,iPar).XLabel.String=parameterName{iPar};
end
for iPar = 1:nPar
    H(iPar).FaceColor = spfitPaper.appliedColor(5);
    for jPar = 1:nPar
        S(iPar,jPar).Color = spfitPaper.appliedColor(5);
        S(iPar, jPar).MarkerSize = 3;
        grid(ax(iPar,jPar),'on');
    end
end
fig.Units               = 'centimeters';
fig.Position(3)         = opts.size;
fig.Position(4)         = opts.size;
set(fig.Children, ...
    'FontName',     opts.font, ...
    'FontSize',     opts.font_size);
set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset+0.02))
fig.PaperPositionMode   = 'auto';
print -dpng -r600 half_hour_scatter_dw
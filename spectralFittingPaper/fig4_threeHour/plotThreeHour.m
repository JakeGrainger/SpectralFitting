%%% script to plot figure 4 from the paper.
clear
%% load fitted data
load('fittedThreeHour')

%% add timings as extra parameter
nPar = length(trueParameter);

%% get options for plotting
boxOptions = spfitPaper.appliedBoxSettings([1,2,5]);
opts = spfitPaper.appliedPlotSettings();

%% plot figure 4
parameterName = {'\alpha','\omega_p','\gamma','r'};
fittingName = {'Least Squares', 'Bartlett LS', 'Debiased Whittle'};
fig = figure();clf;hold on;
ax = spfitPaper.multipleBox(permute(fittedValues(:,:,:),[2,3,1]),...
    fittingName,boxOptions{:});
for jPar = 1:nPar
    ax(jPar).YLabel.String = parameterName{jPar};
    yline(ax(jPar),trueParameter(jPar),'--');
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
print -dpng -r600 three_hour_box

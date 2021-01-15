clear
load('../../table/LargeExperiment.mat')
%%
% get options for plotting
boxOptions = spfitPaper.appliedBoxSettings([1,2,5]);
opts = spfitPaper.appliedPlotSettings();
parameterName = {'\alpha','\omega_p','\gamma','r'};
fittingName = {'Least Squares', 'Bartlett LS', 'Debiased Whittle'};
for iTruePar = 1:nTruePar
    
    trueParameter = parameterChoice{iTruePar};
    fig = figure();clf;hold on;
    ax = spfitPaper.multipleBox(permute(squeeze(fullFits(:,:,:,iTruePar)),[2,3,1]),...
        fittingName,boxOptions{:});
    for jPar = 1:4
        ax(jPar).YLabel.String = parameterName{jPar};
        yline(ax(jPar),trueParameter(jPar),'--');
        grid(ax(jPar),'on');
    end
    sgtitle(sprintf("\\alpha = %0.2f, \\omega_p = %0.2f, \\gamma = %0.2f, r = %0.2f",trueParameter))%(1),trueParameter(2),trueParameter(3),trueParameter(4)))
    fig.Units               = 'centimeters';
    fig.Position(3)         = opts.size;
    fig.Position(4)         = 20;
    set(fig.Children, ...
        'FontName',     opts.font, ...
        'FontSize',     opts.font_size);
    set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
    fig.PaperPositionMode   = 'auto';
    print(sprintf("halfHour_%d", iTruePar), '-dpng', '-r600' )

end

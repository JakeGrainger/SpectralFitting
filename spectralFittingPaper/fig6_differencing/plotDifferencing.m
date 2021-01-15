%%% script to plot figure 2 and 3 from the paper.
clear
%% load fitted data
load('fittedDifferencing')

%% get options for plotting
boxOptions = spfitPaper.appliedBoxSettings([1,2,5]);
boxOptions{6} = repmat(boxOptions{6},2,1);
boxOptions{4} = 1:6;
opts = spfitPaper.appliedPlotSettings();
parameterName = {'\alpha','\omega_p','\gamma','r'};

%% perform plotting
for ii = 1:2
    fittedValues(:,:,1) = lsValues{ii};
    fittedValues(:,:,2) = blsValues{ii};
    fittedValues(:,:,3) = dwValues{ii};
    fittedValues(:,:,4) = lsDiff{ii};
    fittedValues(:,:,5) = blsDiff{ii};
    fittedValues(:,:,6) = dwDiff{ii};
    fig = figure();clf;hold on;
    ax = spfitPaper.multipleBox(permute(fittedValues,[2,3,1]),{'LS','BLS','DW','LSd','BLSd','DWd'},boxOptions{:});
    for jPar = 1:4
        ax(jPar).YLabel.String = parameterName{jPar};
        yline(ax(jPar),trueParameter(jPar),'--');
        grid(ax(jPar),'on');
        if jPar == 1
            set(ax(jPar),'ylim',[0.1,2]);
        elseif jPar == 2
            set(ax(jPar),'ylim',[0.73,0.89]);
        elseif jPar > 2
            set(ax(jPar),'ylim',[0.5,13.5]);
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
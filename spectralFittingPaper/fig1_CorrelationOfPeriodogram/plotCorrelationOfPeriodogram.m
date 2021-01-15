%%% script to generate figure 1, correlation of the periodogram at
%%% different sampling intervals, with and without differencing.
clear
%% get plot options
opts = spfitPaper.appliedPlotSettings();

%% make model and set parameters
sdf = @spectralFitting.genJONSWAPsdf;
theta = [0.7,0.7,3.3,5];
lb = [0,0,1,1];
ub = Inf(4,1);
mackaySpecModel = spectralFitting.SpectralModel(sdf,theta,lb,ub);
N = 4096;

%% loop and produce figures
for Delta = [1,1/4]
    for difference = [0,1]
        fig = figure();clf;
        mackaySpecModel.useDifferencing = difference;
        colormap(cool)
        plotCorrI(mackaySpecModel,N,Delta)
        xlim = ([0,pi/Delta]); ylim = ([0,pi/Delta]);
        if Delta == 1
            xticks((0:1/2:1)*pi)
            xticklabels(["0","\pi/2","\pi"])
            yticks((0:1/2:1)*pi)
            yticklabels(["0","\pi/2","\pi"])
            if difference == 1
                colorbar('location','westoutside','visible','off')
                colorbar('location','eastoutside','Limits',[0,1])
            else
                colorbar('location','westoutside','visible','off')
                colorbar('location','eastoutside','visible','off')
            end
        elseif Delta == 1/4
            xticks((0:2:4)*pi)
            xticklabels(["0","2\pi","4\pi"])
            yticks((0:2:4)*pi)
            yticklabels(["0"," 2\pi"," 4\pi"])
            colorbar('location','westoutside','visible','off')
            colorbar('location','eastoutside','visible','off')
        end
        xlabel("\omega")
        ylabel("\omega")
        grid on; box on
        fig.Units               = 'centimeters';
        fig.Position(3)         = 8.5;
        fig.Position(4)         = 4.5;
        set(fig.Children, ...
            'FontName',     opts.font, ...
            'FontSize',     opts.font_size);
        set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
        fname = sprintf('corr_sampleFreq-%d_diff-%d',1/Delta,difference);
        print(fname, '-dpng', '-r600' )
    end
end
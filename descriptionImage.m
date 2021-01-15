% code to generate image for package description
sdf = @spectralFitting.genJONSWAPsdf;
theta = [0.7;0.7;3.3;4];
lb = [0;0;1;1];
ub = Inf(4,1);
waveSpecModel = spectralFitting.SpectralModel(sdf,theta,lb,ub);
N = 2048;
Delta = 1;
gp1 = waveSpecModel.simulateGaussianProcess(N,Delta);
%%
opts.size = 6;
opts.font = 'Times';
opts.font_size = 8;
opts.inset = 0.04;
fig = figure();clf;
plotdB(spectralFitting.Periodogram(gp1,Delta),'color',[0,0,0,0.5]);hold on
plotdBEI(waveSpecModel,N,Delta);
grid on;
box on;
xticks(0:pi/4:pi)
xticklabels(["0","\pi/4","\pi/2","3\pi/4","\pi"])
ylim([-50,10])
xlim([0,pi])
fig.Units               = 'centimeters';
fig.Position(3)         = opts.size;
fig.Position(4)         = opts.size;
set(fig.Children, ...
    'FontName',     opts.font, ...
    'FontSize',     opts.font_size);
set(gca,'LooseInset',max(get(gca,'TightInset'), opts.inset))
fname = 'toolboxIcon';
print(fname, '-dpng', '-r600' )
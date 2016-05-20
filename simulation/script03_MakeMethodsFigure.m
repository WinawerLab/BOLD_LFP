
%% get all the correlation values 
clear all
sim_nr      = 1;
set_nr      = 1;
load(['./data/NS_simnr' int2str(sim_nr) '_set' int2str(set_nr)],'NS')

%% 
numNeurons2plot     = 5;
numNeuronsColor     = {};
trial2plot          = 1;
t_int               = 50:250;

figure('Position',[0 0 600 300]);
ax = get(gca); % to get colororder to make bars the same color as the neuronal timeseries

subplot(3,6,1),hold on
for k=1:numNeurons2plot
    plot(k*2+NS.data.bb_inputs(t_int,k,trial2plot))
end
axis tight

subplot(3,6,2),hold on
for k=1:numNeurons2plot
    plot(k*max(abs(NS.params.poisson_g))+NS.data.g_inputs(t_int,k,trial2plot))
end
axis tight

subplot(3,6,3),hold on
for k=1:numNeurons2plot
    plot(k*max(abs(NS.params.poisson_a))+NS.data.a_inputs(t_int,k,trial2plot))
end
axis tight

subplot(3,6,4),hold on
signal2plot = NS.data.bb_inputs+NS.data.g_inputs+NS.data.a_inputs;
for k=1:numNeurons2plot
    plot(k*2+signal2plot(t_int,k,trial2plot))
end
axis tight

subplot(3,6,5),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
for k=1:numNeurons2plot
    plot(k*.6+signal2plot(t_int,k))
end
axis tight

subplot(3,6,6),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
for k=1:numNeurons2plot
    barh(k,mean(signal2plot(:,k).^2,1),'FaceColor',ax.ColorOrder(k,:))
end
ylim([0 numNeurons2plot+1])

subplot(3,6,11),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
plot(sum(signal2plot(t_int,:),2),'k')
axis tight

subplot(3,6,12),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
barh(1,sum(mean(signal2plot.^2,1),2),'k')
ylim([0 numNeurons2plot+1])

subplot(3,6,16),hold on
neuron2plot = 5;
ts_for_fft = squeeze(NS.data.ts(:,neuron2plot,NS.params.trials_save_inputs(trial2plot)));
% apply a window
w = window(@hann,1/NS.params.dt);
ts_for_fft = bsxfun(@times, ts_for_fft, w);
% fft power
pxx = abs(fft(ts_for_fft)).^2;
f = ns_get(NS,'f');
f_lims = [find(f>=5,1) find(f>=100,1)];
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'Color',ax.ColorOrder(neuron2plot,:))
set(gca,'XScale','log','YScale','log')
axis tight

subplot(3,6,17),hold on
pxx = ns_get(NS,'power');
f = ns_get(NS,'f');
f_lims = [find(f>=5,1) find(f>=100,1)];
signal2plot = pxx(:,NS.params.trials_save_inputs(trial2plot));
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'k')
set(gca,'XScale','log','YScale','log')
axis tight

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/ns_MethodsFig01'])
% print('-depsc','-r300',['./figures/ns_MethodsFig01'])

%% plot one neuron
%% 

neuron2plot         = 5;
numNeuronsColor     = {};
trial2plot          = 1;
t_int               = 70:240;

figure('Position',[0 0 800 50]);
ax = get(gca); % to get colororder to make bars the same color as the neuronal timeseries

subplot(1,7,1),hold on
plot(1:length(t_int),zeros(size(t_int)),'k')
plot(NS.data.bb_inputs(t_int,neuron2plot,trial2plot))
axis tight

subplot(1,7,2),hold on
plot(1:length(t_int),zeros(size(t_int)),'k')
plot(NS.data.g_inputs(t_int,neuron2plot,trial2plot))
axis tight

subplot(1,7,3),hold on
plot(1:length(t_int),zeros(size(t_int)),'k')
plot(NS.data.a_inputs(t_int,neuron2plot,trial2plot))
axis tight

subplot(1,7,4),hold on
plot(1:length(t_int),zeros(size(t_int)),'k')
signal2plot = NS.data.bb_inputs+NS.data.g_inputs+NS.data.a_inputs;
plot(signal2plot(t_int,neuron2plot,trial2plot))
axis tight

subplot(1,7,5),hold on
plot(1:length(t_int),zeros(size(t_int)),'k')
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
plot(signal2plot(t_int,neuron2plot))
axis tight

subplot(1,7,6),hold on
signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));
barh(0,mean(signal2plot(:,neuron2plot).^2,1),'FaceColor',ax.ColorOrder(1,:))
ylim([-2 2])
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/ns_MethodsFig02_trial' int2str(NS.params.trials_save_inputs(trial2plot))])
print('-depsc','-r300',['./figures/ns_MethodsFig02_trial' int2str(NS.params.trials_save_inputs(trial2plot))])

figure('Position',[0 0 500 60])
f = ns_get(NS,'f');
f_lims = [find(f>=5,1) find(f>=100,1)];

subplot(1,7,1),hold on
ts_for_fft = NS.data.bb_inputs(:,neuron2plot,trial2plot);
pxx = ns_fftpower(ts_for_fft);
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'Color',ax.ColorOrder(1,:))
axis tight

subplot(1,7,2),hold on
ts_for_fft = NS.data.g_inputs(:,neuron2plot,trial2plot);
pxx = ns_fftpower(ts_for_fft);
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'Color',ax.ColorOrder(1,:))
axis tight

subplot(1,7,3),hold on
ts_for_fft = NS.data.a_inputs(:,neuron2plot,trial2plot);
pxx = ns_fftpower(ts_for_fft);
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'Color',ax.ColorOrder(1,:))
axis tight

subplot(1,7,4),hold on
signal2plot = NS.data.bb_inputs+NS.data.g_inputs+NS.data.a_inputs;
ts_for_fft = signal2plot(:,neuron2plot,trial2plot);
pxx = ns_fftpower(ts_for_fft);
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'Color',ax.ColorOrder(1,:))
axis tight

subplot(1,7,5),hold on
hold on
ts_for_fft = squeeze(NS.data.ts(:,neuron2plot,NS.params.trials_save_inputs(trial2plot)));
pxx = ns_fftpower(ts_for_fft);
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'Color',ax.ColorOrder(1,:))
set(gca,'XScale','log','YScale','log')
axis tight
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/ns_MethodsFig02_pxx_trial' int2str(NS.params.trials_save_inputs(trial2plot))])
print('-depsc','-r300',['./figures/ns_MethodsFig02_pxx_trial' int2str(NS.params.trials_save_inputs(trial2plot))])

%% plot a bunch of neurons

neurons2plot         = 1:10;
numNeuronsColor     = {};
trial2plot          = 1;
t_int               = 70:240;

signal2plot = NS.data.ts(:,:,NS.params.trials_save_inputs(trial2plot));

figure('Position',[0 0 150 120])
hold on
for k=1:length(neurons2plot)
    plot(t_int,1/3*k+signal2plot(t_int,neurons2plot(k)),'k')    
end
axis tight
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/ns_MethodsFig02_moreneurons'])
print('-depsc','-r300',['./figures/ns_MethodsFig02_moreneurons'])

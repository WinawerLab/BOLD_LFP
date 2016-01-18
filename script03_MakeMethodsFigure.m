
%% get all the correlation values 
clear all
sim_nr = 1;
poisson_bb_rg_in = [0 .5;0 1];
load(['../neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set1'],'NS')

%% 
numNeurons2plot = 5;
trial2plot = 1;
t_int = 50:250;
figure('Position',[0 0 600 300]);

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
    barh(k,mean(signal2plot(:,k).^2,1))
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

subplot(3,6,17),hold on
pxx = ns_get(NS,'power');
f = ns_get(NS,'f');
f_lims = [find(f>=5,1) find(f>=100,1)];
signal2plot = pxx(:,NS.params.trials_save_inputs(trial2plot));
plot(f(f_lims(1):f_lims(2)),pxx(f_lims(1):f_lims(2)),'k')
set(gca,'XScale','log','YScale','log')
axis tight

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['../neural_sim_output/figures/ns_MethodsFig01'])
print('-depsc','-r300',['../neural_sim_output/figures/ns_MethodsFig01'])

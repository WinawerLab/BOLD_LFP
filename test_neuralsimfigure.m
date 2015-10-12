%% make initial figures SfN poster

[~,high_bb_trial] = max(NS.trial.poisson_rate_bb);
[~,high_gamma_trial] = max(NS.trial.poisson_rate_g);
[~,high_alpha_trial] = max(NS.trial.poisson_rate_a);

num_exam_neur = 5;

figure('Position',[0 0 600 800])
% plot high BB trial
subplot(5,3,1:2),hold on
rand_neurons_b = randperm(ns_get(NS,'num_broadband'),num_exam_neur);
for k=1:num_exam_neur
    plot(k+NS.data.ts(:,rand_neurons_b(k),high_bb_trial),'Color',[.1 .8 .3])
end
xlim([0 1000]),ylim([0 6])
ylabel([int2str(num_exam_neur) ' neurons'])
title([' bb=' num2str(NS.trial.poisson_rate_bb(high_bb_trial)),...
    ' g=' num2str(NS.trial.poisson_rate_g(high_bb_trial)),...
    ' a=' num2str(NS.trial.poisson_rate_a(high_bb_trial))])
set(gca,'YTick',[1:5],'YTickLabel',[])

% plot high G trial
subplot(5,3,4:5),hold on
rand_neurons_g = ns_get(NS,'num_broadband')+randperm(ns_get(NS,'num_gamma'),num_exam_neur);
for k=1:num_exam_neur
    plot(k+NS.data.ts(:,rand_neurons_g(k),high_gamma_trial),'b')
end
xlim([0 1000]),ylim([0 6])
ylabel([int2str(num_exam_neur) ' neurons'])
title([' bb=' num2str(NS.trial.poisson_rate_bb(high_gamma_trial)),...
    ' g=' num2str(NS.trial.poisson_rate_g(high_gamma_trial)),...
    ' a=' num2str(NS.trial.poisson_rate_a(high_gamma_trial))])
set(gca,'YTick',[1:5],'YTickLabel',[])

% plot high alpha trial
subplot(5,3,7:8),hold on
rand_neurons_a = randperm(ns_get(NS,'num_broadband'),num_exam_neur);
for k=1:num_exam_neur
    plot(k+NS.data.ts(:,rand_neurons_a(k),high_alpha_trial),'Color',[.8 .1 .6])
end
xlim([0 1000]),ylim([0 6])
ylabel([int2str(num_exam_neur) ' neurons'])
title([' bb=' num2str(NS.trial.poisson_rate_bb(high_alpha_trial)),...
    ' g=' num2str(NS.trial.poisson_rate_g(high_alpha_trial)),...
    ' a=' num2str(NS.trial.poisson_rate_a(high_alpha_trial))])
set(gca,'YTick',[1:5],'YTickLabel',[])

subplot(5,3,10:11),hold on
plot(140+zeros(1,size(NS.data.ts,1)),'k')
plot(70+zeros(1,size(NS.data.ts,1)),'k')
plot(0+zeros(1,size(NS.data.ts,1)),'k')
plot(140+sum(NS.data.ts(:,:,high_bb_trial),2),'Color',[.1 .8 .3],'LineWidth',2)
plot(70+sum(NS.data.ts(:,:,high_gamma_trial),2),'b','LineWidth',2)
plot(0+sum(NS.data.ts(:,:,high_alpha_trial),2),'Color',[.8 .1 .6],'LineWidth',2)
axis tight
set(gca,'YTick',[])
ylabel('sum across 600 neurons')

subplot(5,3,3),hold on
for k=1:num_exam_neur
    plot(var(NS.data.ts(:,rand_neurons_b(k),high_bb_trial)),k,'*','Color',[.1 .8 .3])
end
xlim([0 .17]),ylim([0 6])
xlabel('variance'),ylabel([int2str(num_exam_neur) ' neurons'])
set(gca,'YTick',[1:5],'YTickLabel',[])

subplot(5,3,6),hold on
for k=1:num_exam_neur
    plot(var(NS.data.ts(:,rand_neurons_g(k),high_gamma_trial)),k,'*','Color',[0 0 1])
end
xlim([0 .17]),ylim([0 6])
xlabel('variance'),ylabel([int2str(num_exam_neur) ' neurons'])
set(gca,'YTick',[1:5],'YTickLabel',[])

subplot(5,3,9),hold on
for k=1:num_exam_neur
    plot(var(NS.data.ts(:,rand_neurons_a(k),high_alpha_trial)),k,'*','Color',[.8 .1 .6])
end
xlim([0 .17]),ylim([0 6])
xlabel('variance'),ylabel([int2str(num_exam_neur) ' neurons'])
set(gca,'YTick',[1:5],'YTickLabel',[])

% subplot(5,3,[12 15]),hold on
% f = ns_get(NS, 'f');
% w = window(@hann,1/NS.params.dt);
% ts_for_fft = squeeze(mean(NS.data.ts(:,:,high_bb_trial),2));
% ts_for_fft = bsxfun(@times, ts_for_fft, w);
% pxx1 = abs(fft(ts_for_fft)).^2;
% plot(f(f<200),pxx1(f<200),'Color',[.1 .8 .3])
% 
% ts_for_fft = squeeze(mean(NS.data.ts(:,:,high_gamma_trial),2));
% ts_for_fft = bsxfun(@times, ts_for_fft, w);
% pxx2 = abs(fft(ts_for_fft)).^2;
% plot(f(f<200),pxx2(f<200),'Color',[0 0 1])
% 
% ts_for_fft = squeeze(mean(NS.data.ts(:,:,high_alpha_trial),2));
% ts_for_fft = bsxfun(@times, ts_for_fft, w);
% pxx3 = abs(fft(ts_for_fft)).^2;
% plot(f(f<200),pxx3(f<200),'Color',[.8 .1 .6])
% 
% set(gca, 'XScale', 'log', 'YScale', 'log')
% xlabel ('Frequency')
% ylabel('Power')
% xlim([1 200]);
% ylim([min([pxx1; pxx2; pxx3]) max([pxx1; pxx2; pxx3])])

subplot(5,3,12),hold on
plot(var(sum(NS.data.ts(:,:,high_bb_trial),2)),140,'*','Color',[.1 .8 .3])
plot(var(sum(NS.data.ts(:,:,high_gamma_trial),2)),70,'*','Color',[0 0 1])
plot(var(sum(NS.data.ts(:,:,high_alpha_trial),2)),0,'*','Color',[.8 .1 .6])
ylim([-30 200])
xlabel('variance')
set(gca,'YTick',[])

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['../figures/simulated_timeseries'])
print('-depsc','-r300',['../figures/simulated_timeseries'])

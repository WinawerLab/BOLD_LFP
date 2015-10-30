
%% ECOG BOLD simulation RUN A LOOP

% Purpose: Simulate neural data - time varying membrane potentials - and
% then ask whether the simulated BOLD signal and various metrics of the
% simulated field potentials are correlated.

%%%%% vary broadband 
sim_nr = 1;
poisson_bb_rg_in = [0 0;0 0.05;0 .1;0 .2;0 .3;0 .4;0 .5;0 1];
sim_nr = 2;
poisson_bb_rg_in = [0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1];
sim_nr = 3;
poisson_bb_rg_in = [0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2];
tic

for k=1:size(poisson_bb_rg_in,1)

    % Set default parameters.
    NS = neural_sim_defaults; %disp(NS.params)
    % Change these to change the simulation. 

    NS = ns_set(NS, 'poisson_bb_rg', poisson_bb_rg_in(k,:));
    NS = ns_set(NS, 'poisson_baseline', .1); 
    
    if sim_nr == 1
        NS = ns_set(NS, 'poisson_g_rg', [0 .5]);
        NS = ns_set(NS, 'poisson_a_rg', [0 .5]); 
    else
        NS = ns_set(NS, 'poisson_g_rg', [0 .5]);
        NS = ns_set(NS, 'poisson_a_rg', [0 .5]); 
    end
    % Assign expected values of broadband and gamma levels for each stimulus class and each trial
    NS = ns_make_trial_struct(NS); %disp(NS.trial)

    % Simulate. This will produce a time series for each neuron in each trial
    NS = ns_simulate_data(NS); 

    % Convert the neural time series into instrument measures
    NS = ns_neural2instruments(NS); %disp(NS.data)

    % Compute the correlations between different instrument measures 
    NS = ns_summary_statistics(NS); %disp(NS.stats)
    
    NS.data.ts = single(NS.data.ts); % to reduce size
    
    save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(k) ],'NS')
    disp(['done simulation ' int2str(k) ' of ' int2str(size(poisson_bb_rg_in,1))])
    toc
end

%% get all the correlation values 
clear all
% sim_nr = 1;
% poisson_bb_rg_in = [0 0;0 0.05;0 .1;0 .2;0 .3;0 .4;0 .5;0 1];
sim_nr = 2;
poisson_bb_rg_in = [0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1];
% sim_nr = 3;
% poisson_bb_rg_in = [0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2;0 .2];

bold_corr=zeros(size(poisson_bb_rg_in,1),5);
bold_beta=zeros(size(poisson_bb_rg_in,1),6);
r_all = zeros(size(poisson_bb_rg_in,1),501);
out.bold_vals = zeros(size(poisson_bb_rg_in,1),8);
out.bb_vals = zeros(size(poisson_bb_rg_in,1),8);
out.gamma_vals = zeros(size(poisson_bb_rg_in,1),8);
out.alpha_vals = zeros(size(poisson_bb_rg_in,1),8);
out.mean_vals = zeros(size(poisson_bb_rg_in,1),8);

for k=1%:size(poisson_bb_rg_in,1)
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(k) ],'NS')
    
    %- to correlate bold and neurophys
    alpha_avg = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'alpha')));
    bb_avg  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bb')));
    gamma_avg  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'gamma')));
    bold_avg  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bold')));
    all_mean_data = squeeze(mean(sum(NS.data.ts(:,:,:),2),1));
    mean_avg = zeros(max(NS.trial.condition_num),1);
    for m=1:max(NS.trial.condition_num)+1
        mean_avg(m) = mean(all_mean_data(NS.trial.condition_num==m-1));
    end
    mean_avg = zscore(mean_avg);
    
    % put in output
    out.bold_vals(k,:) = bold_avg;
    out.bb_vals(k,:) = bb_avg;
    out.gamma_vals(k,:) = gamma_avg;
    out.alpha_vals(k,:) = alpha_avg;
    out.mean_vals(k,:) = mean_avg;
    
    % do stats
    stats = regstats(bold_avg,bb_avg);
    bold_corr(k,1) = stats.rsquare;
    bold_beta(k,1) = stats.beta(2);

    stats = regstats(bold_avg,gamma_avg);
    bold_corr(k,2) = stats.rsquare;
    bold_beta(k,2) = stats.beta(2);
    
    stats = regstats(bold_avg,alpha_avg);
    bold_corr(k,3) = stats.rsquare;
    bold_beta(k,3) = stats.beta(2);
    
    stats = regstats(bold_avg,mean_avg);
    bold_corr(k,4) = stats.rsquare;
    bold_beta(k,4) = stats.beta(2);

    stats = regstats(bold_avg,[bb_avg alpha_avg]);
    bold_corr(k,5) = stats.rsquare;
    bold_beta(k,[5:6]) = stats.beta(2:3);
    
    %- correlate across all frequencies
    bold_avg  = ns_mean_by_stimulus(NS, ns_get(NS, 'bold'));
    ts_for_fft = squeeze(mean(ns_get(NS, 'ts'),2));
    
    % fft settings just like in the other correlation
    fft_w = 250; % window width
    fft_ov = 0; % overlap
    srate = 1/ns_get(NS,'dt');
    % initialize 
    [~,f] = pwelch(ts_for_fft(:,1),fft_w,fft_ov,srate,srate);
    pxx_all = zeros(size(ts_for_fft,2),length(f));
    % loop over trials
    for m = 1:size(ts_for_fft,2)
        pxx_all(m,:) = pwelch(ts_for_fft(:,m),fft_w,fft_ov,srate,srate);
    end
    % mean power by stimulus
    pxx_avg = zeros(NS.params.num_conditions,length(f));
    for m=1:NS.params.num_conditions
        pxx_avg(m,:) = mean(pxx_all(NS.trial.condition_num==m-1,:));
    end
    % power change
    pxx_change = bsxfun(@minus, log10(pxx_avg), log10(pxx_avg(1,:)));
    r = zeros(size(pxx_change,2),1);
    for m=1:size(pxx_change,2)
        r(m) = corr(bold_avg,pxx_change(:,m));
    end
    r_all(k,:) = r;
    
end

%% Plot beta and r2 for each model separately
figure('Position',[0 0 1000 300])
for k=1:size(poisson_bb_rg_in,1)
    subplot(2,size(poisson_bb_rg_in,1),k),hold on
    bar(1-.3,bold_corr(k,1),.3,'r')
    bar(1,bold_corr(k,2),.3,'y')
    bar(1+.3,bold_corr(k,3),.3,'b')
    bar(1+.6,bold_corr(k,5),.3,'FaceColor',[1 0 1])
    xlim([.5 2]),ylim([0 1])
    title(['bb [0 ' int2str(100*poisson_bb_rg_in(k,2)) ']/100' ])
end
xlim([.5 1.9]),ylim([0 1])
for k=1:size(poisson_bb_rg_in,1)
    subplot(2,size(poisson_bb_rg_in,1),size(poisson_bb_rg_in,1)+k),hold on
    bar(1-.3,bold_beta(k,1),.1,'FaceColor',[.9 .9 .9])
    bar(1,bold_beta(k,2),.1,'FaceColor',[.6 .6 .6])
    bar(1+.3,bold_beta(k,3),.1,'FaceColor',[.2 .3 .3])
    bar([1+.6],bold_beta(k,5),.10,'FaceColor',[.9 .9 .9])
    bar([1+.8],bold_beta(k,6),.10,'FaceColor',[.2 .3 .3])
    xlim([.5 2]),ylim([-1.5 1.5])
end

set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/r_varyBB_sim6_k1'])
% print('-depsc','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/r_varyBB_sim6_k1'])


%% Plot average r2 and beta and std across models
figure('Position',[0 0 200 300])    
subplot(2,1,1),hold on
bar(1,mean(bold_corr(:,1)),1,'r')
errorbar(1,mean(bold_corr(:,1)),std(bold_corr(:,1))/sqrt(length(bold_corr(:,1))),'k')
bar(2,mean(bold_corr(:,2)),1,'y')
errorbar(2,mean(bold_corr(:,2)),std(bold_corr(:,2))/sqrt(length(bold_corr(:,2))),'k')
bar(3,mean(bold_corr(:,3)),1,'b')
errorbar(3,mean(bold_corr(:,3)),std(bold_corr(:,3))/sqrt(length(bold_corr(:,3))),'k')
bar(4,mean(bold_corr(:,5)),1,'FaceColor',[1 0 1])
errorbar(4,mean(bold_corr(:,5)),std(bold_corr(:,5))/sqrt(length(bold_corr(:,5))),'k')
xlim([.5 4.5]),ylim([0 1])

subplot(2,1,2),hold on
bar(1-.3,mean(bold_beta(:,1)),.3,'FaceColor',[.9 .9 .9])
errorbar(1-.3,mean(bold_beta(:,1)),std(bold_beta(:,1))/sqrt(length(bold_beta(:,1))),'k')
bar(2.3-.3,mean(bold_beta(:,2)),.3,'FaceColor',[.6 .6 .6])
errorbar(2.3-.3,mean(bold_beta(:,2)),std(bold_beta(:,2))/sqrt(length(bold_beta(:,2))),'k')
bar(3.6-.3,mean(bold_beta(:,3)),.3,'FaceColor',[.2 .3 .3])
errorbar(3.6-.3,mean(bold_beta(:,3)),std(bold_beta(:,3))/sqrt(length(bold_beta(:,3))),'k')
bar(4-.3,mean(bold_beta(:,5)),.3,'FaceColor',[.9 .9 .9])
errorbar(4-.3,mean(bold_beta(:,5)),std(bold_beta(:,5))/sqrt(length(bold_beta(:,5))),'k')
bar(4.6-.3,mean(bold_beta(:,6)),.3,'FaceColor',[.2 .3 .3])
errorbar(4.6-.3,mean(bold_beta(:,6)),std(bold_beta(:,6))/sqrt(length(bold_beta(:,6))),'k')
xlim([.5 4.5])
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/r_varyBB_sim' int2str(sim_nr) '_av'])
print('-depsc','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/r_varyBB_sim' int2str(sim_nr) '_av'])

%%
figure('Position',[0 0 200 150]),hold on
plot(f,zeros(size(f)),'Color',[.5 .5 .5])
plot(f,r_all(4,:),'Color',[.5 .5 .5],'LineWidth',1)
plot(f,mean(r_all),'k','LineWidth',2)
xlabel('Frequency'),ylabel('correlation with BOLD (r)')
xlim([0 200])
ylim([-1 1])

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/allfreq_corr1val_nr' int2str(sim_nr)])
print('-depsc','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/allfreq_corr1val_nr' int2str(sim_nr)])

%% load a set and plot stuff 
%%
%% plot inputs table
%%

clear all

% sim_nr = 5;
% poisson_bb_rg_in = [0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1];
% s=1;
sim_nr = 6;
poisson_bb_rg_in = [0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5];
s=6; % S4 cond 1-5; S6 c1-6; S8 c1-3;

bold_corr=zeros(size(poisson_bb_rg_in,1),5);
bold_beta=zeros(size(poisson_bb_rg_in,1),6);
r_all = zeros(size(poisson_bb_rg_in,1),501);

load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(s) ],'NS')

% to plot the mean with alpha
all_mean_data = squeeze(mean(sum(NS.data.ts(:,:,:),2),1));
mean_avg = zeros(max(NS.trial.condition_num),1);
for m=1:max(NS.trial.condition_num)+1
    mean_avg(m) = mean(all_mean_data(NS.trial.condition_num==m-1));
end


%%

plot_colors = [0 0 0; jet(NS.params.num_conditions)];

figure('Position',[0 0 400 400])
subplot(7,1,1),hold on
for k=1:length(NS.params.poisson_bb)
    bar(k,NS.params.poisson_baseline+NS.params.poisson_bb(k),'FaceColor',plot_colors(k,:))
end
ylim([.99 1.5])
set(gca,'YTick',[1 1.5])
ylabel('bb level')
subplot(7,1,2),hold on
for k=1:length(NS.params.poisson_g)
    bar(k,0+NS.params.poisson_g_rg(2),'FaceColor',plot_colors(k,:))
end
ylabel('gamma level')
subplot(7,1,3),hold on
for k=1:length(NS.params.poisson_a)
    bar(k,0+NS.params.poisson_a_rg(2),'FaceColor',plot_colors(k,:))
end
ylabel('alpha level')

subplot(7,1,4),hold on
for k=1:length(NS.params.poisson_bb)
    bar(k,0,'FaceColor',plot_colors(k,:))
end
ylabel('bb coh')
subplot(7,1,5),hold on
for k=1:length(NS.params.poisson_g)
    bar(k,NS.params.poisson_g(k)/max(NS.params.poisson_g),'FaceColor',plot_colors(k,:))
end
ylabel('gamma coh')
subplot(7,1,6),hold on
for k=1:length(NS.params.poisson_a)
    bar(k,NS.params.poisson_a(k)/max(NS.params.poisson_a),'FaceColor',plot_colors(k,:))
end
ylabel('alpha coh')

subplot(7,1,7),hold on
for k=1:length(NS.params.poisson_a)
    bar(k,mean_avg(k),.5,'FaceColor',plot_colors(k,:))
end
ylabel('mean')
% bar([1:8],1./(1+NS.params.poisson_a),.5,'c')
% ylabel('mean inputs')

for k=1:7
    subplot(7,1,k)
    box off
    xlim([0 9])
    set(gca,'XTick',[1:8])
end
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/example_inputs' int2str(sim_nr) '_' int2str(s)])
print('-depsc','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/example_inputs' int2str(sim_nr) '_' int2str(s)])

%% figure for SfN poster
%% PLOT

% sim_nr = 1;
% poisson_bb_rg_in = [0 .5; 0 .4;0 .3;0 .2;0 .1;0 0];
% load('/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_varyBB0-0p50_nr1.mat')
% sim_nr = 5;
% poisson_bb_rg_in = [0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1;0 .1];
% k=1;
% load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(k) ],'NS')
sim_nr = 6;
poisson_bb_rg_in = [0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5;0 .5];
s=6;
load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(s) ],'NS')

%%
[bb_avg,bb_std]         = ns_mean_by_stimulus(NS, ns_get(NS, 'bb'));
[bold_avg,bold_std]     = ns_mean_by_stimulus(NS, ns_get(NS, 'bold'));
[lfp_avg,lfp_std]       = ns_mean_by_stimulus(NS, ns_get(NS, 'lfp'));
[gamma_avg,gamma_std]   = ns_mean_by_stimulus(NS, ns_get(NS, 'gamma'));
[alpha_avg,alpha_std]   = ns_mean_by_stimulus(NS, ns_get(NS, 'alpha'));
num_conditions          = ns_get(NS, 'num_conditions');
freq_bb                 = ns_get(NS, 'freq_bb');

fH=figure('Position',[0 0 900 700]); clf; set(fH, 'Color', 'w')
fs = [18 12]; % fontsize

% ---- Plot Spectra for different stimuli -----
subplot(3,3,[1 4]), set(gca, 'FontSize', fs(1));
plot_colors = [0 0 0; jet(num_conditions)];
set(gca, 'ColorOrder', plot_colors); hold all
plot(ns_get(NS, 'f'), ns_mean_by_stimulus(NS, ns_get(NS, 'power')), '-', ...
    freq_bb, exp(ns_get(NS, 'power_law')), 'k-', 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log','XTick',[10 100])
xlabel ('Frequency')
ylabel('Power')
% xlim([min(freq_bb) max(freq_bb)]);
xlim([5 max(freq_bb)]);

% ---- Plot BOLD v ECoG measures ----------------
num_subplots = 4; % broadband; total LFP; gamma; alpha
x_data = {bb_avg, lfp_avg, gamma_avg, alpha_avg};
x_err = {bb_std, lfp_std, gamma_std, alpha_std};
for ii=1:length(x_err)
    x_err{ii} = x_err{ii}./sqrt(NS.params.num_averages);
end
bold_err = bold_std./sqrt(NS.params.num_averages);
xl     = {'broadband', 'Total LFP power', 'Gamma', 'Alpha'};
for ii = 1:num_subplots
    this_subplot = 3*(ii-1)+2;
    subplot(num_subplots,3,this_subplot), set(gca, 'FontSize', fs(2)); hold on
    p = polyfit(x_data{ii}, bold_avg,1);
    error_x = [x_data{ii}-x_err{ii} x_data{ii}+x_err{ii}];
    error_y = [bold_avg-bold_err bold_avg+bold_err];
    plot([x_data{ii} x_data{ii}]',error_y','Color',[.5 .5 .5]);
    plot(error_x',[bold_avg bold_avg]','Color',[.5 .5 .5]);
    scatter(x_data{ii}, bold_avg,40,plot_colors(2:end,:)), axis tight square
    
    hold on; plot(x_data{ii}, polyval(p, x_data{ii}), 'k-', 'LineWidth', 1)
    xlabel(xl{ii}), ylabel('BOLD')
    title(sprintf('r = %4.2f', corr(x_data{ii}, bold_avg)));
end


% ---- Plot BOLD and ECoG measures as function of simulation inputs -----
num_subplots = 3; % broadband; total LFP; gamma; alpha
x_data_name = {'poisson_bb', 'poisson_g', 'poisson_a'};
xl     = {'broadband', 'Gamma', 'Alpha'};

for ii = 1:num_subplots
    this_subplot = 3 * (ii-1)+3;
    subplot(num_subplots,3,this_subplot), set(gca, 'FontSize', fs(2));
    x_data = ns_get(NS, x_data_name{ii});
    
    [~, inds] = sort(x_data);
    plot(...
        x_data(inds), zscore(bold_avg(inds)), 'o-',...
        x_data(inds), zscore(bb_avg(inds)),  'd-',...
        x_data(inds), zscore(lfp_avg(inds)),  's-',...
        x_data(inds), zscore(gamma_avg(inds)), 'x-',...
        x_data(inds), zscore(alpha_avg(inds)), '*-',...
        'LineWidth', 3)
    xlabel(sprintf('%s (inputs)', xl{ii})), ylabel('response (z-scores)')
    if ii == 1
        legend({'BOLD', 'Broadband',  'LFP power',  'Gamma', 'Alpha'},...
            'Location', 'Best', 'Box', 'off')
    end
end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/powerspectra_example_set' int2str(sim_nr) '_' int2str(s)])
% print('-depsc','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/powerspectra_example_set' int2str(sim_nr) '_' int2str(s)])
%%
figure('Position',[0 0 300 300]),hold on
for k=1:8
    bar(k,bold_avg(k),'FaceColor',plot_colors(k,:))
end
set(gca,'XTick',[1:8])
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/powerspectra_example_BOLDset' int2str(sim_nr) '_' int2str(s)])
print('-depsc','-r300',['/Volumes/DoraBigDrive/github/neural_sim_output/figures/powerspectra_example_BOLDset' int2str(sim_nr) '_' int2str(s)])

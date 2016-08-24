
%% ECOG BOLD simulation two main simulations

% Purpose: Simulate neural data - time varying membrane potentials - and
% then ask whether the simulated BOLD signal and various metrics of the
% simulated field potentials are correlated.

%%%% TODO: check alpha - influence by bb (set bb to zero/variable)

sim_nr = 3;

% make conditions of input:
nr_conds = 8;
out = script_make_uncorrelated_conditions(nr_conds);

tic
clear ns_params
for k=4%1:size(out.poisson_bb,1) % number of settings

    % Set default parameters
    NS = neural_sim_defaults; %disp(NS.params)

    NS = ns_set(NS, 'save_inputs', 1);
    NS = ns_set(NS, 'num_neurons', 200); 
    NS = ns_set(NS, 'poisson_baseline', .2); 
    NS = ns_set(NS, 'poisson_bb',out.poisson_bb(k,:));
    NS = ns_set(NS, 'poisson_g',out.poisson_g(k,:));
    NS = ns_set(NS, 'poisson_a',out.poisson_a(k,:));
    NS = ns_set(NS, 'coherence_bb',out.coherence_bb(k,:));
    NS = ns_set(NS, 'coherence_g',out.coherence_g(k,:));
    NS = ns_set(NS, 'coherence_a',out.coherence_a(k,:));
    
    % set bb, gamma and alpha values for all trials
    NS = ns_make_trial_struct(NS);
    
    % if save_inputs choose trials to save, data can get big if saving a lot
    NS = ns_set(NS, 'trials_save_inputs',[1 length(NS.trial.poisson_bb)]);

    % Simulate. This will produce a time series for each neuron in each trial
    NS = ns_simulate_data(NS); 

    % Convert the neural time series into instrument (LFP/BOLD) measures
    NS = ns_neural2instruments(NS); %disp(NS.data)

    % Get spectra and broadband, gamma, alpha and mean from the LFP
    NS = ns_analyse_lfp(NS); %disp(NS.data)
    
    % Compute the correlations between different instrument measures 
    NS = ns_summary_statistics(NS); %disp(NS.stats)
        
    NS.data.ts = single(NS.data.ts); % to reduce size
    
%     save(['../data/NS_simnr' int2str(sim_nr) '_set' int2str(k) ],'NS')
    save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(k) ],'NS')
    ns_params{k} = NS.params;
    disp(['done simulation ' int2str(k) ' of ' int2str(size(out.poisson_bb,1))])
    toc
end

%%
%% figures for checks:

clear all

sim_nr = 3; % simulation number (general settings change)

for param_set=4%1:3%1:10 % parameter set (bb, gamma, alpha)

load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(param_set) ],'NS')

bold_avg        = median(NS.data.bold_bs,2);
bb_avg          = median(NS.data.bb,2);
gamma_avg       = median(NS.data.gamma,2);
alpha_avg       = median(NS.data.alpha,2);
bold_ci         = [quantile(NS.data.bold_bs,.16,2) quantile(NS.data.bold_bs,.84,2)];
bb_ci           = [quantile(NS.data.bb,.16,2) quantile(NS.data.bb,.84,2)];
gamma_ci        = [quantile(NS.data.gamma,.16,2) quantile(NS.data.gamma,.84,2)];
alpha_ci        = [quantile(NS.data.alpha,.16,2) quantile(NS.data.alpha,.84,2)];

num_conditions  = ns_get(NS, 'num_conditions');
f               = ns_get(NS, 'f');

% INPUTS AND OUTPUTS
figure('Position',[0 0 1100 750])  

% ---- Plot Spectra for different stimuli -----
subplot(4,6,1), set(gca, 'FontSize', 10);
plot_colors = [0 0 0; jet(num_conditions-1)];
set(gca, 'ColorOrder', plot_colors); hold all

for k = 1:num_conditions
    plot(f,mean(NS.data.lfp_spectra(:,NS.trial.condition_num==k-1),2),...
        'Color',plot_colors(k,:),'LineWidth',2)
end
set(gca, 'XScale', 'log', 'YScale', 'log','XTick',[10 100],'YTick',[10.^-2 10.^0 10.^2])
xlabel ('Frequency'), ylabel('Power')
xlim([5 max(f)]); %ylim([10.^-3 10.^2]);

% ---- Plot BOLD for different stimuli -----
subplot(4,6,2), set(gca, 'FontSize', 10),hold on
for ii = 1:num_conditions
    bar(ii,bold_avg(ii),'FaceColor',plot_colors(ii,:))
end
ylabel('BOLD')

% ---- Plot BOLD v ECoG measures ----------------
x_data = {bb_avg, gamma_avg, alpha_avg};
x_err = {bb_ci, gamma_ci, alpha_ci};
xl     = {'broadband', 'gamma', 'alpha'};
for ii = 1:length(x_data)
    subplot(4,6,6+ii), hold on
    p = polyfit(x_data{ii}, bold_avg,1);
    error_x = x_err{ii};
    error_y = bold_ci;
    plot([x_data{ii} x_data{ii}]',error_y','-','Color',[.5 .5 .5]);
    plot(error_x',[bold_avg bold_avg]','-','Color',[.5 .5 .5]);
    scatter(x_data{ii}, bold_avg,40,plot_colors), axis tight square
    hold on; plot(x_data{ii}, polyval(p, x_data{ii}), 'k-', 'LineWidth', 1)
    xlabel(xl{ii},'FontSize',10), ylabel('BOLD','FontSize',10)
    title(sprintf('r = %4.2f', corr(x_data{ii}, bold_avg)));
end

% ---- Plot BOLD and ECoG measures as function of simulation inputs -----
x_input     = {'poisson_bb', 'poisson_g', 'poisson_a'};
x_output    = {'bb_avg', 'gamma_avg', 'alpha_avg'};
for ii = 1:length(x_data)
    subplot(4,6,12+ii), hold on
    x = getfield(NS.params,x_input{ii});
    y = eval(x_output{ii});
    for k = 1:length(x)
        plot(x(k),y(k),'.','MarkerSize',20,'Color',plot_colors(k,:))
    end
    xlabel(x_input{ii}), ylabel(x_output{ii})
end

x_input     = {'coherence_bb', 'coherence_g', 'coherence_a'};
x_output    = {'bb_avg', 'gamma_avg', 'alpha_avg'};
for ii = 1:length(x_data)
    subplot(4,6,18+ii), hold on
    x = getfield(NS.params,x_input{ii});
    y = eval(x_output{ii});
    for k = 1:length(x)
        plot(x(k),y(k),'.','MarkerSize',20,'Color',plot_colors(k,:))
    end

    xlabel(x_input{ii}), ylabel(x_output{ii})
end

% ---- Plot BOLD correlation with LFP across frequencies -----
subplot(4,3,3), set(gca, 'FontSize', 10),hold on
r_boot=zeros(size(NS.data.lfp_spectra_bs,2),length(f));
for bs=1:size(NS.data.lfp_spectra_bs,2) % number of bootstraps
    fmri_d=NS.data.bold_bs(:,bs);
    ecog_d=squeeze(NS.data.lfp_spectra_bs(:,bs,:));
    r_boot(bs,:)=corr(fmri_d,ecog_d);
end
plot(f,zeros(size(r_boot)),'k','LineWidth',1)
plot(f,median(r_boot,1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])

% REGRESSION MODEL:
bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};

subplot(4,8,14:16),hold on
for k=1:length(NS.stats)
    cross_val_r2 = median(NS.stats(k).stats(:,3),1);
    bar(k,cross_val_r2,.9,'FaceColor',bar_colors{k})
end
xlim([0 8]),ylim([0 1])
set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb g','a','bb a','g a','bb g a'})

labels_beta={{'bb','',''},{'','g',''},{'bb','g',''},{'','','a'},...
    {'bb','','a'},{'','g','a'},{'bb','g','a'}};
labels_index={[1],[2],[1 2],[3],[1 3],[2 3],[1 2 3]};
bb_g_a_color={[.9 .9 .9],[.6 .6 .6],[.3 .3 .3]};

% get minmax beta for plotting max y
min_y = 0;
max_y = 0;
for k=1:length(NS.stats)
    min_y = min([min_y median(NS.stats(k).beta(:,2:end),1)]);
    max_y = max([max_y median(NS.stats(k).beta(:,2:end),1)]);
end

for k=1:length(NS.stats)
    xl_ind=labels_index{k};
    subplot(4,3*length(NS.stats),2*3*length(NS.stats)+2*length(NS.stats)+k),hold on 

    temp_beta=median(NS.stats(k).beta(:,2:end),1);
    for m=1:length(xl_ind)
        bar(xl_ind(m),temp_beta(m),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
    end

%     xlim([.5 3.5]),ylim([min_y-.1 max_y+.1])
    xlim([.5 3.5]),ylim([-1 1])
    set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-1:.2:1.5],'YTickLabel',[])
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/set' int2str(param_set) '_inputs_outputs'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/set' int2str(param_set) '_inputs_outputs'])

% close all

end

%% figure of the mean and std of the ts

name_x_in = 'poisson_a';
% name_x_in = 'coherence_g';

eval(['x_in = NS.params.' name_x_in;]);

figure
mean_ts = zeros(num_conditions,NS.params.num_neurons);
for k = 1:num_conditions
    % mean over time and condition
    mean_ts(k,:) = mean(mean(NS.data.ts(:,:,NS.trial.condition_num==k-1),1),3);
end

subplot(2,1,1)
plot(x_in,mean(mean_ts,2),'k.','MarkerSize',20)
[r,p]=corr(x_in',mean(mean_ts,2));
xlabel(['input ' name_x_in])
ylabel('timeseries mean')
title(['r.^2=' num2str(r.^2)])

subplot(2,1,2)
plot(x_in,std(mean_ts,[],2),'k.','MarkerSize',20)
[r,p]=corr(x_in',std(mean_ts,[],2));
xlabel(['input ' name_x_in])
ylabel('timeseries std')
title(['r.^2=' num2str(r.^2)])





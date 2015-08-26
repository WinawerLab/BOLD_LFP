%% ECOG BOLD simulation

% Purpose: Simulate neural data - time varying membrane potentials - and
% then ask whether the simulated BOLD signal and various metrics of the
% simulated field potentials are correlated.

% Set default parameters.
NS = neural_sim_defaults; disp(NS.params)
% Change these to change the simulation. 
% For example:  NS = ns_set(NS, 'num_neurons', 300);
% Or, NS = ns_set(NS, 'num_neurons', 50); NS = ns_set(NS, 'num_experiments', 2);
% NS = ns_set(NS, 'num_averages', 10);
NS = ns_set(NS, 'poisson_bb_rg', [0 .3]);
NS = ns_set(NS, 'poisson_g_rg', [0 .3]);
NS = ns_set(NS, 'gamma_coh', 1);

% Assign expected values of broadband and gamma levels for each stimulus class and each trial
NS = ns_make_trial_struct(NS); disp(NS.trial)

% Simulate. This will produce a time series for each neuron in each trial
NS = ns_simulate_data(NS); 

% Convert the neural time series into instrument measures
NS = ns_neural2instruments(NS); disp(NS.data)

% Compute the correlations between different instrument measures 
NS = ns_summary_statistics(NS); disp(NS.stats)

%% PLOT
bb_avg    = ns_mean_by_stimulus(NS, ns_get(NS, 'bb'));
bold_avg  = ns_mean_by_stimulus(NS, ns_get(NS, 'bold'));
lfp_avg   = ns_mean_by_stimulus(NS, ns_get(NS, 'lfp'));
gamma_avg = ns_mean_by_stimulus(NS, ns_get(NS, 'gamma'));
alpha_avg = ns_mean_by_stimulus(NS, ns_get(NS, 'alpha'));
num_conditions = ns_get(NS, 'num_conditions');
freq_bb = ns_get(NS, 'freq_bb');

fH=figure(300); clf; set(fH, 'Color', 'w')
fs = [18 12]; % fontsize

% ---- Plot Spectra for different stimuli -----
subplot(1,3,1), set(gca, 'FontSize', fs(1));
plot_colors = [0 0 0; jet(num_conditions)];
set(gca, 'ColorOrder', plot_colors); hold all
plot(ns_get(NS, 'f'), ns_mean_by_stimulus(NS, ns_get(NS, 'power')), '-', ...
    freq_bb, exp(ns_get(NS, 'power_law')), 'k-', 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel ('Frequency')
ylabel('Power')
xlim([min(freq_bb) max(freq_bb)]);

% ---- Plot BOLD v ECoG measures ----------------
num_subplots = 4; % broadband; total LFP; gamma; alpha
x_data = {bb_avg, lfp_avg, gamma_avg, alpha_avg};
xl     = {'broadband', 'Total LFP power', 'Gamma', 'Alpha'};
for ii = 1:num_subplots
    this_subplot = 3*(ii-1)+2;
    subplot(num_subplots,3,this_subplot), set(gca, 'FontSize', fs(2)); hold on
    p = polyfit(x_data{ii}, bold_avg,1);
    scatter(x_data{ii}, bold_avg), axis tight square
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

return
%%


% Summarize R2 values acorss experiments and measurement types (BOLD v
% broadband, total LFP, or gamma)
if num_experiments > 1
    figure(101);clf
    subplot(3,1,1)
    hist(R2_bb_bold, -.05:.1:1.05); hold on
    plot([1 1]*median(R2_bb_bold), get(gca, 'YLim'), 'k-')
    title(sprintf('Broadband v BOLD R^2, %4.2f',median(R2_bb_bold))), xlim([0 1])
    
    subplot(3,1,2)
    hist(R2_lfp_bold, -.05:.1:1.05); hold on
    plot([1 1]*median(R2_lfp_bold), get(gca, 'YLim'), 'k-')
    title(sprintf('LFP v BOLD R^2, %4.2f',median(R2_lfp_bold))), xlim([0 1])
    
    subplot(3,1,3)
    hist(R2_gamma_bold, -.05:.1:1.05); hold on
    plot([1 1]*median(R2_gamma_bold), get(gca, 'YLim'), 'k-')
    title(sprintf('gamma v BOLD R^2, %4.2f',median(R2_gamma_bold))), xlim([0 1])
end
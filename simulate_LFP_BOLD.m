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

fH=figure(300); clf;

% ---- Plot Spectra for different stimuli -----
subplot(1,3,1)
plot_colors = [0 0 0; jet(num_conditions)];
set(gca, 'ColorOrder', plot_colors); hold all
plot(ns_get(NS, 'f'), ns_mean_by_stimulus(NS, ns_get(NS, 'power')), '-', ...
    freq_bb, exp(ns_get(NS, 'power_law')), 'k-', 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel ('Frequency')
ylabel('Power')
xlim([min(freq_bb) max(freq_bb)]);

% ---- Plot BOLD v ECoG measures ----------------
subplot(3,3,2)
p = polyfit(bb_avg, bold_avg,1);
scatter(bb_avg, bold_avg), axis tight square
hold on; plot(bb_avg, polyval(p, bb_avg), 'k-', 'LineWidth', 1)
xlabel('broadband'), ylabel('BOLD')
title(sprintf('r = %4.2f', corr(bb_avg, bold_avg)));

subplot(3,3,5)
p = polyfit(lfp_avg, bold_avg,1);
scatter(lfp_avg, bold_avg), axis tight square
hold on; plot(lfp_avg, polyval(p, lfp_avg), 'k-', 'LineWidth', 1)
xlabel('Total LFP power'), ylabel('BOLD')
title(sprintf('r = %4.2f', corr(lfp_avg(:), bold_avg)));

subplot(3,3,8)
p = polyfit(gamma_avg, bold_avg,1);
scatter(gamma_avg, bold_avg), axis tight square
hold on; plot(gamma_avg, polyval(p, gamma_avg), 'k-', 'LineWidth', 1)
xlabel('Gamma'), ylabel('BOLD')
title(sprintf('r = %4.2f', corr(gamma_avg, bold_avg)));

% % sanity check: plot (BOLD + COV) v LFP: they should be identical
% figure(5), clf
% scatter(lfp, bold+cov)
% xlabel('LFP'), ylabel('BOLD + COV')
% xl = get(gca, 'XLim');
% set(gca, 'YLim', xl)
% hold on;
% plot(xl, xl, 'k--'); hold off


% ---- Plot BOLD and ECoG measures as function of simulation inputs -----
subplot(2,3,3)
poisson_bb = ns_get(NS, 'poisson_bb');
poisson_g  = ns_get(NS, 'poisson_g');

[~, inds] = sort(poisson_bb);
plot(...
    poisson_bb(inds), zscore(bold_avg(inds)), 'o-',...
    poisson_bb(inds), zscore(bb_avg(inds)),  'd-',...
    poisson_bb(inds), zscore(lfp_avg(inds)),  's-',...
    poisson_bb(inds), zscore(gamma_avg(inds)), 'x-', 'LineWidth', 3)
xlabel('bb level (inputs)'), ylabel('response (z-scores)')
legend({'BOLD', 'Broadband',  'LFP power',  'Gammma'}, 'Location', 'Best', 'Box', 'off')

subplot(2,3,6)
[~, inds] = sort(poisson_g);
plot(...
    poisson_g(inds), zscore(bold_avg(inds)), 'o-',...
    poisson_g(inds), zscore(bb_avg(inds)),  'd-',...
    poisson_g(inds), zscore(lfp_avg(inds)),  's-',...
    poisson_g(inds), zscore(gamma_avg(inds)), 'x-', 'LineWidth', 3)
xlabel('gamma level (inputs)'), ylabel('response (z-scores)')
%legend({'BOLD', 'Broadband',  'LFP power',  'Gammma'}, 'Location', 'Best', 'Box', 'off')

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
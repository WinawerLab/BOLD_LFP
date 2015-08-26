function NS = ns_summary_statistics(NS)
%  Summarize data from simulations
bb      = ns_get(NS, 'bb');
bold    = ns_get(NS, 'bold');
lfp     = ns_get(NS, 'lfp');
gamma   = ns_get(NS, 'gamma');
alpha   = ns_get(NS, 'alpha');

conditions = ns_get(NS, 'condition_num');
% Average trials across repeats of the same stimuli
for ii = 1:ns_get(NS, 'num_conditions')
    inds = conditions == ii-1;
    lfp_avg(ii, :)   = mean(lfp(inds,:));
    bold_avg(ii, :)  = mean(bold(inds,:));
    bb_avg(ii, :)    = mean(bb(inds,:));
    gamma_avg(ii, :) = mean(gamma(inds,:)) - mean(bb(inds,:));
    alpha_avg(ii, :) = mean(alpha(inds,:)) - mean(bb(inds,:));        
end
 

% Compute R2 values for each LFP measure against BOLD, for each experiment
R2_bb_bold    = diag(corr(bb_avg, bold_avg)).^2;
R2_lfp_bold   = diag(corr(lfp_avg, bold_avg)).^2;
R2_gamma_bold = diag(corr(gamma_avg, bold_avg)).^2;
R2_alpha_bold = diag(corr(alpha_avg, bold_avg)).^2;

NS = ns_set(NS, 'R2', R2_bb_bold, 'bb_bold');
NS = ns_set(NS, 'R2', R2_lfp_bold, 'lfp_bold');
NS = ns_set(NS, 'R2', R2_gamma_bold, 'gamma_bold');
NS = ns_set(NS, 'R2', R2_alpha_bold, 'alpha_bold');

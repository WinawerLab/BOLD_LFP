function NS = ns_make_trial_struct(NS,fixed_bb,fixed_g,fixed_a,fixed_g_coh)
% Set the trial-by-trial expected values for broadband and gamma values
% NS = ns_make_trial_struct(NS)

% Broadband and gamma levels are both condition/stimulus dependent. The
% minimum and maximum level (defined in units of Poisson rates) are already
% specified as simulation parameters. Here we linearly space the levels
% across the unique conditions/stimuli
num_conditions = ns_get(NS, 'num_conditions');
poisson_bb_rg  = ns_get(NS, 'poisson_bb_rg');
poisson_g_val  = ns_get(NS, 'poisson_g_val');
poisson_a_rg  = ns_get(NS, 'poisson_a_rg');
coherence_g_rg = ns_get(NS, 'gamma_coh_rg');

% get values random or pre-specified:
if ~exist('fixed_bb', 'var') || isempty(fixed_bb) % broadband
    poisson_bb      = [0 linspace(poisson_bb_rg(1),poisson_bb_rg(2), num_conditions-1)];
else
    poisson_bb      = fixed_bb;     
end
if ~exist('fixed_g', 'var') || isempty(fixed_g) % gamma amplitude
    poisson_g      = poisson_g_val*ones(size(poisson_bb));
else
    poisson_g      = fixed_g;     
end
if ~exist('fixed_g_coh', 'var') || isempty(fixed_g_coh) % gamma coherence
    coherence_g    = (linspace(coherence_g_rg(1), coherence_g_rg(2),  num_conditions-1)).^2;
else
    coherence_g    = fixed_g_coh;
end
if ~exist('fixed_a', 'var') || isempty(fixed_a) % alpha amplitude
    poisson_a      = (linspace(poisson_a_rg(1), poisson_a_rg(2),  num_conditions-1)).^2;
else
    poisson_a      = fixed_a;
end

% Decorrelate expected broadband and gamma levels across conditions /
% stimuli. We would like the measures to be decorrelated across trials 
% because if one of the measures correlates better with BOLD than the
% other, it will be easiest to detect this if the two measures are not
% themselves correlated. We decorrelate empirically: 
% Suppose the broadband levels for 4 stimuli are 1, 2, 4, and 8. 
% And suppose the gamma levels are also 1, 2, 4, 8 (but not necessarily in
% this order). We decorrelate by repeatedly scrambling the order of
% stimulus to level and choosing the order that is most orthogonal. There
% is probably a smarter way to do this with linear algebra. But anyway,
% this method is fast, so who cares.

if ~exist('fixed_g_coh', 'var') || isempty(fixed_g_coh) % only decorrelate if no fixed input
    b = zeros(100, num_conditions);
    c = zeros(num_conditions,1);
    for ii = 1:100;
        tmp =  [0 coherence_g(randperm(num_conditions-1))];
        b(ii,:) = tmp;
        c(ii) = corr(poisson_bb', tmp');
    end
    [~, ind]  = min(abs(c));
    coherence_g = b(ind,:);
end
% Decorrelate alpha (in same way as bb/gamma)
% Randomize alpha levels across trials
if ~exist('fixed_a', 'var') || isempty(fixed_a) % only decorrelate if no fixed input
    b = zeros(1000, num_conditions);
    c = zeros(num_conditions,2); % correlation with bb and gamma
    for ii = 1:1000;
        tmp =  [max(poisson_a) poisson_a(randperm(num_conditions-1))];
        b(ii,:) = tmp;
        c(ii,1) = corr(poisson_bb', tmp');
        c(ii,2) = corr(coherence_g', tmp');
    end
    [~, ind]  = min(sum(abs(c),2));
    poisson_a = b(ind,:);
end

% At this point we have one gamma level and one broadband level for each
% unique condition / stimulus. We now assign levels to all trials (which
% can include repeated conditions / stimuli).
num_trials = ns_get(NS, 'num_trials');
vectorize = @(data,repeats) reshape(repmat(data, repeats, 1), num_trials, 1);
NS = ns_set(NS, 'condition_num',   vectorize(0:num_conditions-1, ns_get(NS, 'num_averages')));

NS = ns_set(NS, 'poisson_rate_bb', vectorize(poisson_bb, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'poisson_rate_g', vectorize(poisson_g, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'poisson_rate_a', vectorize(poisson_a, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'coherence_rate_g',  vectorize(coherence_g, ns_get(NS, 'num_averages')));

NS = ns_set(NS, 'poisson_bb', poisson_bb);
NS = ns_set(NS, 'poisson_g', poisson_g);
NS = ns_set(NS, 'poisson_a', poisson_a);
NS = ns_set(NS, 'coherence_g', coherence_g);

% Visualize
% figure;  
% plot(NS.trial.poisson_rate_bb, 'x');  hold on;
% plot(NS.trial.poisson_rate_g, 'd');
% plot(NS.trial.poisson_rate_a, '*');
% legend('Broadband poisson rate', 'Gamma poisson rate', 'Alpha poisson rate')
% set(gca, 'XTick', (0:num_conditions)*ns_get(NS, 'num_averages'), 'XGrid', 'on'); 
return





end
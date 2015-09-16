function NS = ns_make_trial_struct(NS)
% Set the trial-by-trial expected values for broadband and gamma values
% NS = ns_make_trial_struct(NS)

% Broadband and gamma levels are both condition/stimulus dependent. The
% minimum and maximum level (defined in units of Poisson rates) are already
% specified as simulation parameters. Here we linearly space the levels
% across the unique conditions/stimuli
num_conditions = ns_get(NS, 'num_conditions');
poisson_bb_rg  = ns_get(NS, 'poisson_bb_rg');
poisson_g_rg   = ns_get(NS, 'poisson_g_rg', 'scaled');
poisson_a_rg   = ns_get(NS, 'poisson_a_rg', 'scaled');
poisson_bb     = [0 linspace(poisson_bb_rg(1),poisson_bb_rg(2), num_conditions-1)];
poisson_g      = linspace(poisson_g_rg(1), poisson_g_rg(2),  num_conditions-1);
poisson_a      = linspace(poisson_a_rg(1), poisson_a_rg(2),  num_conditions-1);

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

b = zeros(100, num_conditions);
c = zeros(num_conditions,1);
for ii = 1:100;
    tmp =  [0 poisson_g(randperm(num_conditions-1))];
    b(ii,:) = tmp;
    c(ii) = corr(poisson_bb', tmp');
end
[~, ind]  = min(abs(c));
poisson_g = b(ind,:);

% Decorrelate alpha (in same way as bb/gamma
% Randomize alpha levels across trials
b = zeros(100, num_conditions);
c = zeros(num_conditions,2); % correlation with bb and gamma
for ii = 1:100;
    tmp =  [0 poisson_a(randperm(num_conditions-1))];
    b(ii,:) = tmp;
    c(ii,1) = corr(poisson_bb', tmp');
    c(ii,2) = corr(poisson_g', tmp');
end
[~, ind]  = min(sum(abs(c),2));
poisson_a = b(ind,:);

% At this point we have one gamma level and one broadband level for each
% unique condition / stimulus. We now assign levels to all trials (which
% can include repeated conditions / stimuli).
num_trials = ns_get(NS, 'num_trials');
vectorize = @(data,repeats) reshape(repmat(data, repeats, 1), num_trials, 1);
NS = ns_set(NS, 'condition_num',   vectorize(0:num_conditions-1, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'poisson_rate_bb', vectorize(poisson_bb, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'poisson_rate_g',  vectorize(poisson_g, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'poisson_rate_a',  vectorize(poisson_a, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'poisson_g', poisson_g);
NS = ns_set(NS, 'poisson_bb', poisson_bb);
NS = ns_set(NS, 'poisson_a', poisson_a);

% Visualize
figure;  plot(NS.trial.condition_num, 'o'); hold on;
plot(NS.trial.poisson_rate_bb, 'x'); 
plot(NS.trial.poisson_rate_g, 'd');
plot(NS.trial.poisson_rate_a, '*');
legend('Condition number', 'Broadband poisson rate', 'Gamma poisson rate', 'Alpha poisson rate')

return





end
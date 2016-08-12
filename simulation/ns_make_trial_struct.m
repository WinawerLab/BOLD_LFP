function NS = ns_make_trial_struct(NS)
% Set the trial-by-trial expected values for broadband and gamma values
% NS = ns_make_trial_struct(NS)

% At this point we have one bb, gamma and alpha levels and coherences for
% each unique condition / stimulus. We now assign levels to all trials
% (which can include repeated conditions / stimuli).
num_trials      = ns_get(NS, 'num_trials');
num_conditions  = ns_get(NS, 'num_conditions');
poisson_bb      = ns_get(NS, 'poisson_bb');
poisson_g       = ns_get(NS, 'poisson_g');
poisson_a       = ns_get(NS, 'poisson_a');
coherence_bb    = ns_get(NS, 'coherence_bb');
coherence_g     = ns_get(NS, 'coherence_g');
coherence_a     = ns_get(NS, 'coherence_a');

% make function to repeat across conditions:
vectorize = @(data,repeats) reshape(repmat(data, repeats, 1), num_trials, 1);

% repeat the stuff across conditions:
NS = ns_set(NS, 'condition_num',   vectorize(0:num_conditions-1, ns_get(NS, 'num_averages')));

NS = ns_set(NS, 'trial_poisson_bb', vectorize(poisson_bb, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'trial_poisson_g', vectorize(poisson_g, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'trial_poisson_a', vectorize(poisson_a, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'trial_coherence_bb',  vectorize(coherence_bb, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'trial_coherence_g',  vectorize(coherence_g, ns_get(NS, 'num_averages')));
NS = ns_set(NS, 'trial_coherence_a',  vectorize(coherence_a, ns_get(NS, 'num_averages')));



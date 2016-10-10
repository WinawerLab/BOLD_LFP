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

% make sure that these are all rows:
if size(poisson_bb,1)>size(poisson_bb,2)
   poisson_bb = poisson_bb'; 
end
if size(poisson_g,1)>size(poisson_g,2)
   poisson_g = poisson_g'; 
end
if size(poisson_a,1)>size(poisson_a,2)
   poisson_a = poisson_a'; 
end
if size(coherence_bb,1)>size(coherence_bb,2)
   coherence_bb = coherence_bb'; 
end
if size(coherence_g,1)>size(coherence_g,2)
   coherence_g = coherence_g'; 
end
if size(coherence_a,1)>size(coherence_a,2)
   coherence_a = coherence_a'; 
end

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



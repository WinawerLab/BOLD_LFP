function NS = ns_neural2instruments_alternativeBOLD(NS)
% 
% NS = ns_neural2instruments(NS)
%
% calculate the LFP and BOLD from the simulated membrane potentials

num_trials      = ns_get(NS, 'num_trials');
all_ts          = ns_get(NS, 'ts');
num_experiments = ns_get(NS, 'num_experiments');

% We assume the BOLD signal is the sum of the energy demand on each neuron,
% and that the energy demand of a neuron is proportional to the variance in
% its membrane potential. We assume the LFP is the sum of the time varying
% membrane potential in each neuron. We then take the variance in the LFP
% as the metric of the LFP response. We also compute the sum of the
% covariance, excluding the diagonal of the covariance matrix
bold_fun  = @(x) sum(mean(abs(x))); %  bold_fun  = @(x) sum(var(x));
lfp_fun   = @(x) sum(x,2);

% intialize variables to summarize trial data
lfp   = zeros(size(all_ts,1),num_trials,num_experiments);
bold  = zeros(num_trials,num_experiments);

for sim_number = 1:num_experiments
    
    ts = all_ts(:,:,:,sim_number);
        
    % compute each summary metric for each trial
    for ii = 1:num_trials
        lfp(:,ii,sim_number) = lfp_fun(ts(:,:,ii));
        bold(ii,sim_number)  = bold_fun(ts(:,:,ii));
    end

end

NS = ns_set(NS, 'lfp', lfp);
NS = ns_set(NS, 'bold', bold);



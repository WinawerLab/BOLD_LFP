function NS = ns_neural2instruments(NS)
%
% NS = ns_neural2instruments(NS)

freq_bb         = ns_get(NS, 'freq_bb');
freq_gamma      = ns_get(NS, 'freq_gamma');
num_trials      = ns_get(NS, 'num_trials');
all_ts          = ns_get(NS, 'ts');
num_experiments = ns_get(NS, 'num_experiments');
data_power      = ns_get(NS, 'power'); 
power_law       = ns_get(NS, 'power_law');

% We assume the BOLD signal is the sum of the energy demand on each neuron,
% and that the energy demand of a neuron is proportional to the variance in
% its membrae potential. We assume the LFP is the sum of the time varying
% membrane potential in each neuron. We then take the variance in the LFP
% as the metric of the LFP response. We also compute the sum of the
% covariance, excluding the diagonal of the covariance matrix
bold_fun  = @(x) sum(mean(x.^2)); % sum(var(x));
lfp_fun   = @(x) var(sum(x,2));
bb_fun    = @(x, power_law_baseline) exp(mean(log(x) - power_law_baseline));
gamma_fun = @(x) exp(mean(log(x)));

% intialize variables to summarize trial data
lfp   = zeros(num_trials,num_experiments);
bold  = zeros(num_trials,num_experiments);
bb    = zeros(num_trials,num_experiments);
gamma = zeros(num_trials,num_experiments);


for sim_number = 1:num_experiments
    
    ts = all_ts(:,:,:,sim_number);
            
    % compute each summary metric for each trial
    for ii = 1:num_trials
        lfp(ii,sim_number)   = lfp_fun(ts(:,:,ii));
        bold(ii,sim_number)  = bold_fun(ts(:,:,ii));
        bb(ii,sim_number)    = bb_fun(data_power(freq_bb,ii), power_law(:,sim_number));
        gamma(ii,sim_number) = gamma_fun(data_power(freq_gamma,ii));
    end
        
end

NS = ns_set(NS, 'lfp', lfp);
NS = ns_set(NS, 'bold', bold);
NS = ns_set(NS, 'bb', bb);
NS = ns_set(NS, 'gamma', gamma);



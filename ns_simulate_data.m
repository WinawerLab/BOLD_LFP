function NS = ns_simulate_data(NS)

% Get variables from NS struct before simulating
t                = ns_get(NS, 't');
num_trials       = ns_get(NS, 'num_trials');
num_broadband    = ns_get(NS, 'num_broadband');
num_gamma        = ns_get(NS, 'num_gamma');
poisson_rate_bb  = ns_get(NS, 'poisson_rate_bb');
poisson_baseline = ns_get(NS, 'poisson_baseline');
poisson_rate_g   = ns_get(NS, 'poisson_rate_g');
gamma_filter     = ns_get(NS, 'gamma_filter');

% Initialize the time series array, which will hold data for all neurons at
% all time points in all trials across all experiments
ts  = zeros(length(t), ns_get(NS, 'num_neurons'), num_trials, ...
    ns_get(NS, 'num_experiments'));

% Loop over experiments
for sim_number = 1:ns_get(NS, 'num_experiments')
    
    fprintf('\n'); drawnow();
    
    % generate the simulated data: time x neuron x trial
    switch ns_get(NS, 'simulate_method')
        
        case 'FAST' % a random walk
            % Initialize variables. These variables will store one time
            % series per tiral per neuron in the broadband pool (ts_bb) and
            % one time series per trial per neuron in the gamma pool (ts_g)
            ts_bb = zeros(length(t), num_broadband, num_trials);
            ts_g  = zeros(length(t), num_gamma, num_trials);
            
            drawdots = round((1:10)/10*num_trials);
            for ii = 1:num_trials
                if ismember(ii,drawdots), fprintf('.'); drawnow(); end
                % Generate seed time series with random numbers scaled by
                % the expected broadband level for the given stimulus class
                this_rate       = poisson_rate_bb(ii) + poisson_baseline;
                ts_bb(:,:,ii)   = randn(length(t), num_broadband)*this_rate;
                
                % We do the same for the gamma time series subject to the
                % additional constraint that the gamma signal is correlated
                % across neurons with coherence <gamma_coh>. Note that we
                % generate a time series with power at all frequencies for
                % the gamma signal, and then band pass filter the result.
                % This is more efficient than generate a band limited
                % random walk.
                mu        = zeros(1,num_gamma);
                sigma     = eye(num_gamma) + (1-eye(num_gamma))* ns_get(NS, 'gamma_coh');
                pre_gamma = mvnrnd(mu,sigma,length(t));
                
                
                gamma_signal    = poisson_rate_g(ii)*filter(gamma_filter, pre_gamma);
                baseline        = randn(length(t), num_gamma)*poisson_baseline;
                ts_g(:,:,ii)    = bsxfun(@plus, gamma_signal, baseline);
                
            end
            
            % The time series of random numbers are turned into random
            % walks by computing the cummualtive sum. This computation
            % turns a white noise stimulus into a brown noise stimulus.
            ts_bb = cumsum(ts_bb);
            ts_g  = cumsum(ts_g);
            % Combine the gamma and broadband populations into a single
            % neural pool.
            ts(:,:,:,sim_number) = cat(2, ts_bb, ts_g);            
            
        case 'SLOW' % a more complicated leaky integrator
            error('Slow method not yet implemented')
    end
    
    NS = ns_set(NS, 'ts', ts);
    
    
    
end

return


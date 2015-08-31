function NS = ns_simulate_data(NS)

% Get variables from NS struct before simulating
t                = ns_get(NS, 't');
num_trials       = ns_get(NS, 'num_trials');
num_broadband    = ns_get(NS, 'num_broadband');
num_gamma        = ns_get(NS, 'num_gamma');
num_alpha        = ns_get(NS, 'num_alpha');
poisson_rate_bb  = ns_get(NS, 'poisson_rate_bb');
poisson_baseline = ns_get(NS, 'poisson_baseline');
poisson_rate_g   = ns_get(NS, 'poisson_rate_g');
poisson_rate_a   = ns_get(NS, 'poisson_rate_a');
gamma_filter     = ns_get(NS, 'gamma_filter');
alpha_filter     = ns_get(NS, 'alpha_filter');

% Initialize the time series array, which will hold data for all neurons at
% all time points in all trials across all experiments
ts  = zeros(length(t), ns_get(NS, 'num_neurons'), num_trials, ...
    ns_get(NS, 'num_experiments'));

% Loop over experiments
for sim_number = 1:ns_get(NS, 'num_experiments')
    
    fprintf('\n'); drawnow();
    
    % Initialize variables. These variables will store one time
    % series per tiral per neuron in the broadband pool (ts_bb) and
    % one time series per trial per neuron in the gamma pool (ts_g)
    ts_bb = zeros(length(t), num_broadband, num_trials);
    ts_g  = zeros(length(t), num_gamma, num_trials);
    
    fprintf('[%s]: Simulating time series ', mfilename);
    drawdots = round((1:10)/10*num_trials);
    
    % generate the simulated data: time x neuron x trial
    switch ns_get(NS, 'simulate_method')
        case 'FAST' % a random walk
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
                % This is more efficient than generating a band limited
                % random walk.
                mu        = zeros(1,num_gamma);
                sigma     = eye(num_gamma) + (1-eye(num_gamma))* ns_get(NS, 'gamma_coh');
                gamma_inputs = mvnrnd(mu,sigma,length(t));
                
                
                gamma_signal    = poisson_rate_g(ii)*filter(gamma_filter, gamma_inputs);
                baseline        = randn(length(t), num_gamma)*poisson_baseline;
                ts_g(:,:,ii)    = bsxfun(@plus, gamma_signal, baseline);
                
            end
            
            fprintf('Done\n');
            % The time series of random numbers are turned into random
            % walks by computing the cummualtive sum. This computation
            % turns a white noise stimulus into a brown noise stimulus.
            ts_bb = cumsum(ts_bb);
            ts_g  = cumsum(ts_g);
            
        case 'SLOW' % a more complicated leaky integrator
            alpha = 0.010;
            dt    = ns_get(NS, 'dt');
            
            % generate data for one trial at a time
            for ii = 1:num_trials
                if ismember(ii,drawdots), fprintf('.'); drawnow(); end
                
                %%%%% Broadband inputs
                this_rate = poisson_rate_bb(ii) + poisson_baseline;
                bb_inputs    = randn(length(t), num_broadband)*this_rate;
                
                %%%%% Gamma inputs
                mu           = zeros(1,num_gamma);
                sigma        = eye(num_gamma) + (1-eye(num_gamma))* ns_get(NS, 'gamma_coh');
                gamma_inputs = mvnrnd(mu,sigma,length(t));
                
                gamma_signal    = poisson_rate_g(ii)*filter(gamma_filter, gamma_inputs);
                baseline        = randn(length(t), num_gamma)*poisson_baseline;
                gamma_inputs    = bsxfun(@plus, gamma_signal, baseline);
                
                %%%%% Alpha inputs
                mu           = zeros(1,num_broadband); % if you add offset here it would get filtered out
                sigma        = eye(num_broadband) + (1-eye(num_broadband))* ns_get(NS, 'alpha_coh');
                alpha_inputs = mvnrnd(mu,sigma,length(t));
                
                % get the alpha signal
                alpha_signal = ns_alpha_signal(alpha_inputs,poisson_rate_a(ii),dt,0);
                                
                % do we need to add baseline here? broadband already has baseline
%                 baseline        = randn(length(t), num_broadband)*poisson_baseline;
%                 alpha_inputs    = bsxfun(@plus, alpha_signal, baseline);
                alpha_inputs    = alpha_signal;
                
                % combine broadband and alpha
                bb_inputs = bb_inputs + alpha_inputs;
                
                for jj = 1:length(t)-1
                    
                    % rate of change in current
                    dIdt = (bb_inputs(jj,:) - ts_bb(jj,:,ii)) / alpha;
                    
                    % stepwise change in current
                    dI = dIdt * dt;
                    
                    % current at next time point
                    ts_bb(jj+1,:,ii) = ts_bb(jj,:,ii) + dI;
                                        
                    % rate of change in current
                    dIdt = (gamma_inputs(jj,:) - ts_g(jj,:,ii)) / alpha;
                    
                    % stepwise change in current
                    dI = dIdt * dt;
                    
                    % current at next time point
                    ts_g(jj+1,:,ii) = ts_g(jj,:,ii) + dI;
                end
            end
            fprintf('Done\n');
    end
    
    
    % Combine the gamma and broadband populations into a single
    % neural pool.
    ts(:,:,:,sim_number) = cat(2, ts_bb, ts_g);
    
    NS = ns_set(NS, 'ts', ts);
    
    
    
end

return


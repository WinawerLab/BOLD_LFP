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

% currently not used, I am not sure whether this is necessary
%%%% DESIGN LOWPASS FILTER FOR THE ENVELOPE
band = [3];
srate = 1/NS.params.dt;
low_Rp = 3; low_Rs = 60; % order Butterworth
low_high_p = band(1)*2/srate;
low_high_s = (band(1)+20)*2/srate;
[low_n_band, low_wn_band] = buttord(low_high_p, low_high_s, low_Rp, low_Rs);
[low_bf_b, low_bf_a] = butter(low_n_band, low_wn_band,'low');

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
    
    tau = 0.010;              % time constant of leaky integrator (seconds)
    dt    = ns_get(NS, 'dt'); % step size for simulation, in seconds
    
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
        
        gamma_inputs    = poisson_rate_g(ii)*filter(gamma_filter, gamma_inputs);
        baseline        = randn(length(t), num_gamma)*poisson_baseline;
        gamma_inputs    = bsxfun(@plus, gamma_inputs, baseline);
        
        %%%%% Alpha inputs
        mu              = zeros(1,num_broadband); % if you add offset here it would get filtered out
        sigma           = eye(num_broadband) + (1-eye(num_broadband))* ns_get(NS, 'alpha_coh');
        alpha_inputs    = mvnrnd(mu,sigma,length(t));
        alpha_inputs    = poisson_rate_a(ii)*filter(alpha_filter, alpha_inputs);        
        alpha_envelope  = abs(hilbert(alpha_inputs));
        % I am not sure that this is necessary:
%         alpha_envelope  = filtfilt(low_bf_b, low_bf_a, alpha_envelope); % lowpass filter the envelope
        alpha_inputs    = alpha_inputs + alpha_envelope;

        % combine broadband and alpha
        bb_inputs = bb_inputs + alpha_inputs;
                
        % Leaky integrator loop 
        for jj = 1:length(t)-1
            
            % rate of change in current
            dIdt = (bb_inputs(jj,:) - ts_bb(jj,:,ii)) / tau;
             
            % stepwise change in current
            dI = dIdt * dt;
            
            % current at next time point
            ts_bb(jj+1,:,ii) = ts_bb(jj,:,ii) + dI;
            
            % rate of change in current
            dIdt = (gamma_inputs(jj,:) - ts_g(jj,:,ii)) / tau;
            
            % stepwise change in current
            dI = dIdt * dt;
            
            % current at next time point
            ts_g(jj+1,:,ii) = ts_g(jj,:,ii) + dI;
        end
        
    ts_bb(:,:,ii) = ts_bb(:,:,ii) + alpha_envelope;
        
    end
    fprintf('Done\n');        
    
    
    % Combine all populations into a single neural pool.
    ts(:,:,:,sim_number) = cat(2, ts_bb, ts_g);
    
%     % subtract baseline mean
     baseline_ts = ts(:,:,ns_get(NS,'baseline_trials'),sim_number);
     ts = ts - mean(baseline_ts(:));
    
    % Store the time series in the NS struct
    NS = ns_set(NS, 'ts', ts);
    
    
    
end

return


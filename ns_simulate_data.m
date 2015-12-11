function NS = ns_simulate_data(NS)

% Get variables from NS struct before simulating
t                = ns_get(NS, 't');
num_trials       = ns_get(NS, 'num_trials');
num_neurons      = ns_get(NS, 'num_neurons');
poisson_baseline = ns_get(NS, 'poisson_baseline');
poisson_rate_bb  = ns_get(NS, 'poisson_rate_bb');
poisson_rate_g   = ns_get(NS, 'poisson_rate_g');
poisson_rate_a   = ns_get(NS, 'poisson_rate_a');
coherence_rate_g = ns_get(NS, 'coherence_rate_g'); % call just coherence
coherence_rate_a = ns_get(NS, 'coherence_rate_a');
gamma_filter     = ns_get(NS, 'gamma_filter');
alpha_filter     = ns_get(NS, 'alpha_filter');
lowpass_filter   = ns_get(NS, 'lowpass_filter');
length_zero_pad  = length(t);

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
    ts_bb = zeros(length(t), num_neurons, num_trials);
%     ts_g  = zeros(length(t), num_neurons, num_trials);
    
    fprintf('[%s]: Simulating time series ', mfilename);
    drawdots = round((1:10)/10*num_trials);
    
    % generate the simulated data: time x neuron x trial
    
    tau = 0.010;              % time constant of leaky integrator (seconds)
    dt    = ns_get(NS, 'dt'); % step size for simulation, in seconds
    
    % generate data for one trial at a time
    for ii = 1:num_trials
        if ismember(ii,drawdots), fprintf('.'); drawnow(); end
        
        %%%%% Broadband inputs
        this_rate       = poisson_rate_bb(ii) + poisson_baseline;
        bb_inputs       = randn(length(t), num_neurons)*this_rate;
        
        %%%%% Gamma inputs
        mu              = zeros(1,num_neurons);
        sigma           = eye(num_neurons) + (1-eye(num_neurons))* coherence_rate_g(ii);
        gamma_inputs    = mvnrnd(mu,sigma,length(t));
        gamma_inputs    = padarray(gamma_inputs, [length_zero_pad 0], 0, 'both'); % zero pad 
        gamma_inputs    = poisson_rate_g(ii)*filtfilt(gamma_filter, gamma_inputs);
        gamma_inputs    = gamma_inputs(length(t)+1:2*length(t),:); % remove zero pad 
%         baseline        = randn(length(t), num_neurons)*poisson_baseline;
%         gamma_inputs    = bsxfun(@plus, gamma_inputs, baseline);
        
        %%%%% Alpha inputs with pulses
        t_total = dt:dt:num_neurons;
        f = (0:length(t_total)-1)/max(t_total);
        alpha_pulses = zeros(size(t_total));
        next_pulse = 0; 
        while next_pulse < max(t_total)
            next_pulse = next_pulse + 0.1+randn*.020;
            [~, idx] = min(abs(next_pulse-t_total));
            alpha_pulses(idx) = .1*rand+1;
        end
        alpha_pulses = alpha_pulses*coherence_rate_a(ii);
        t_h = dt:dt:.300;
%         h = exp(-(t_h-.075).^2/(2*.015^2));
        h = exp(-(t_h-.075).^2/(2*.02^2));
        alpha_signal = -conv(alpha_pulses, h, 'same');
%         alpha_signal = alpha_signal(1:length(t_total));
%         alpha_signal = filtfilt(lowpass_filter, alpha_signal); % lowpass to reduce harmonics
        alpha_signal_epoched = reshape(alpha_signal, length(t), num_neurons);
        alpha_inputs = alpha_signal_epoched;
        
%         mu              = zeros(1,num_neurons); % if you add offset here it would get filtered out
%         sigma           = eye(num_neurons) + (1-eye(num_neurons))*0;%* coherence_rate_a(ii);
%         alpha_inputs    = mvnrnd(mu,sigma,length(t));
%         alpha_inputs    = padarray(alpha_inputs, [length_zero_pad 0], 0, 'both'); % zero pad 
% %         alpha_inputs    = poisson_rate_a(ii)*filtfilt(alpha_filter, alpha_inputs);       
%         alpha_inputs    = coherence_rate_a(ii)*filtfilt(alpha_filter, alpha_inputs);       
%         alpha_envelope  = abs(hilbert(alpha_inputs)); % computed here on every neuron
%         alpha_envelope  = alpha_envelope(length(t)+1:2*length(t),:); % remove zero pad 
%         alpha_inputs    = alpha_inputs(length(t)+1:2*length(t),:); % remove zero pad 
% %         alpha_inputs    = alpha_inputs + (2/(1+coherence_rate_a(ii)))*alpha_envelope;
%         % make alpha into negative pulses
%         alpha_inputs    = -(alpha_envelope+alpha_inputs);
                
        % combine broadband and alpha
        bb_inputs = bb_inputs + alpha_inputs + gamma_inputs;
                
        % Leaky integrator loop 
        for jj = 1:length(t)-1
            
            % rate of change in current
            dIdt = (bb_inputs(jj,:) - ts_bb(jj,:,ii)) / tau;
             
            % stepwise change in current
            dI = dIdt * dt;
            
            % current at next time point
            ts_bb(jj+1,:,ii) = ts_bb(jj,:,ii) + dI;
            
%             % rate of change in current
%             dIdt = (gamma_inputs(jj,:) - ts_g(jj,:,ii)) / tau;
%             
%             % stepwise change in current
%             dI = dIdt * dt;
%             
%             % current at next time point
%             ts_g(jj+1,:,ii) = ts_g(jj,:,ii) + dI;
        end
                
    end
    fprintf('Done\n');        
    
    
    % Combine all populations into a single neural pool.
%     ts(:,:,:,sim_number) = cat(2, ts_bb, ts_g);
    ts(:,:,:,sim_number) = ts_bb;
    
    % subtract baseline mean
    baseline_ts = ts(:,:,ns_get(NS,'baseline_trials'),sim_number);
    ts = ts - mean(baseline_ts(:));
    
    % Store the time series in the NS struct
    NS = ns_set(NS, 'ts', ts);
    
    
    
end

return


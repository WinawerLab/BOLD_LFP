function NS = ns_simulate_data(NS)

% Get variables from NS struct before simulating
t                       = ns_get(NS, 't');
num_trials              = ns_get(NS, 'num_trials');
num_neurons             = ns_get(NS, 'num_neurons');
poisson_baseline        = ns_get(NS, 'poisson_baseline');
poisson_bb              = ns_get(NS, 'trial_poisson_bb');
poisson_g               = ns_get(NS, 'trial_poisson_g');
poisson_a               = ns_get(NS, 'trial_poisson_a');
coherence_bb            = ns_get(NS, 'trial_coherence_bb'); 
coherence_g             = ns_get(NS, 'trial_coherence_g'); 
coherence_a             = ns_get(NS, 'trial_coherence_a');
gamma_filter            = ns_get(NS, 'gamma_filter');
alpha_filter            = ns_get(NS, 'alpha_filter');
length_zero_pad         = length(t);

% Initialize the time series array, which will hold data for all neurons at
% all time points in all trials across all experiments
ts  = zeros(length(t), ns_get(NS, 'num_neurons'), num_trials, ...
    ns_get(NS, 'num_experiments'));

% Initialize the alpha covariance calculated from the alpha inputs:
NS.trial.coherence_a_data = NaN(num_trials,1);

% check whether some inputs need to be saved, and which trials
save_inputs = ns_get(NS, 'save_inputs');
save_inputs_trials = ns_get(NS, 'trials_save_inputs');

% initiate variables to save
save_bb_inputs = zeros(length(t), ns_get(NS, 'num_neurons'),length(save_inputs_trials),ns_get(NS, 'num_experiments'));
save_g_inputs = zeros(length(t), ns_get(NS, 'num_neurons'),length(save_inputs_trials),ns_get(NS, 'num_experiments'));
save_a_inputs = zeros(length(t), ns_get(NS, 'num_neurons'),length(save_inputs_trials),ns_get(NS, 'num_experiments'));

% Loop over experiments
for sim_number = 1:ns_get(NS, 'num_experiments')
    counter_save_trials = 0;
    fprintf('\n'); drawnow();
    
    % Initialize variables. These variables will store one time
    % series per trial per neuron  
    ts_integrated = zeros(length(t), num_neurons, num_trials);
    
    fprintf('[%s]: Simulating time series ', mfilename);
    drawdots = round((1:10)/10*num_trials);
    
    % parameters
    tau         = 0.010;            % time constant of leaky integrator (seconds)
    dt          = ns_get(NS, 'dt'); % step size for simulation, in seconds
    input_dc    = .25;              % DC added, this needs to be the size of the maximum alpha
    if 2*input_dc<max(NS.trial.poisson_a)
        fprintf(['WARNING: do not make alpha sent mean signal below zero \n' ...
        ' change maximum poisson alpha range to <' num2str(input_dc*2)  '\n'])
    end
    
    % generate the simulated data: time x neuron x trial one trial at a time
    for ii = 1:num_trials
        if ismember(ii,drawdots), fprintf('.'); drawnow(); end
        
        %%%%% Broadband inputs
        this_rate       = poisson_bb(ii) + poisson_baseline;
        mu              = input_dc + zeros(1,num_neurons);
        sigma           = eye(num_neurons) + (1-eye(num_neurons)) * coherence_bb(ii);
        bb_inputs       = input_dc + this_rate * (mvnrnd(mu,sigma,length(t))-input_dc);
        
        %%%%% Gamma inputs
        mu              = zeros(1,num_neurons);
        sigma           = eye(num_neurons) + (1-eye(num_neurons)) * coherence_g(ii);
        gamma_inputs    = mvnrnd(mu,sigma,length(t));
        % filter gamma for the gamma frequencies
        gamma_inputs    = padarray(gamma_inputs, [length_zero_pad 0], 0, 'both'); % zero pad 
        gamma_inputs    = poisson_g(ii)*filtfilt(gamma_filter, gamma_inputs);
        gamma_inputs    = gamma_inputs(length(t)+1:2*length(t),:); % remove zero pad 
        
        %%%%% Alpha inputs
        mu              = zeros(1,num_neurons);
        sigma           = eye(num_neurons) + (1-eye(num_neurons)) * coherence_a(ii);
        alpha_inputs    = mvnrnd(mu,sigma,length(t));
        % filter alpha for the alpha frequencies
        alpha_inputs    = padarray(alpha_inputs, [length_zero_pad 0], 0, 'both'); % zero pad 
        alpha_inputs    = filtfilt(alpha_filter, alpha_inputs);
        % now add the hilbert envelope to set the mean to zero
        alpha_envelope  = abs(hilbert(alpha_inputs)); % computed here on every neuron
        alpha_envelope  = alpha_envelope(length(t)+1:2*length(t),:); % remove zero pad 
        alpha_inputs    = alpha_inputs(length(t)+1:2*length(t),:); % remove zero pad         
        alpha_inputs    = alpha_inputs + alpha_envelope;
        alpha_inputs    = -poisson_a(ii) * alpha_inputs;
              
        % combine broadband, gamma and alpha
        summed_inputs = bb_inputs + gamma_inputs + alpha_inputs;
          
        % write out some of the raw inputs if saving
        if save_inputs == 1
            if ismember(ii,save_inputs_trials)
                counter_save_trials = counter_save_trials+1;
                save_bb_inputs(:,:,counter_save_trials,sim_number) = single(bb_inputs);
                save_g_inputs(:,:,counter_save_trials,sim_number) = single(gamma_inputs);
                save_a_inputs(:,:,counter_save_trials,sim_number) = single(alpha_inputs);
            end
        end
        
        % Leaky integrator loop 
        ts_integrated(1,:,ii) = summed_inputs(1,:);
        for jj = 1:length(t)-1
            % rate of change in current
            dIdt = (summed_inputs(jj,:) - ts_integrated(jj,:,ii)) / tau;
            % stepwise change in current
            dI = dIdt * dt;
            % current at next time point
            ts_integrated(jj+1,:,ii) = ts_integrated(jj,:,ii) + dI;
        end
        
    end
    fprintf('Done\n');        
    
    % Combine all populations into a single neural pool.
    ts(:,:,:,sim_number) = ts_integrated;
    
    % Store the time series in the NS struct
    NS = ns_set(NS, 'ts', ts);

    % Store the inputs to the structure if required
    if save_inputs == 1
        NS = ns_set(NS, 'bb_input', save_bb_inputs);
        NS = ns_set(NS, 'g_input', save_g_inputs);
        NS = ns_set(NS, 'a_input', save_a_inputs);
    end
end

return


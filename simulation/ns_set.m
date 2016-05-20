function NS = ns_set(NS, param, val, varargin)
% Set parameters for neural simulation
% NS = ns_set(NS, param, val, varargin)

% Big switch
switch lower(param)
    case 'simulate_method'
        % Two ways to simulate the neural time series:
        %   'SLOW': leaky integrator simulation
        %   'FAST': random walk
        NS.params.simulate_method = val;
        
    case 'save_inputs'
        % save alpha, gamma and broadband inputs before leaky integrator
        NS.params.save_inputs = val;

    case 'trials_save_inputs'
        % save alpha, gamma and broadband inputs before leaky integrator
        NS.params.trials_save_inputs = val;
        
    case 'num_neurons'
        % The number of neurons influences how much the LFP power can deviate from
        % BOLD. This is because the sum of random signals grows with sqrt(n), and
        % synchronous signals grow with n. The larger the value of n, the more the
        % LFP signal (but not the BOLD signal!) will be dominated by the
        % synchronous signals.
        NS.params.num_neurons = val;
        
    case 'num_conditions'
        % Number of stimulus conditions. Each condition is defined by a
        % broadband and a gamma amplitude.
        NS.params.num_conditions = val;
        
    case 'num_averages'
        % Number of trials with the same condition in each experiment
        NS.params.num_averages =  val;
        
    case 'num_experiments'
        % Number of experiments. (> 1 only if we want to compute a distribution of
        % R2 values)
        NS.params.num_experiments = val;
        
    case 'dt'
        % One second per trial, ms sampling
        NS.params.dt =  val;
        
    case 'trial_length'
        NS.params.trial_length = val;
        
    case 'gamma_range'
        % Band limits of spectrum that define gamma band response (Hz)
        NS.params.gamma_range = val;
        
    case 'alpha_range'
        % Band limits of spectrum that define alpha band response (Hz)
        NS.params.alpha_range =  val;
        
    case 'poisson_baseline'
        % Baseline Poisson rate (arbitrary units)
        NS.params.poisson_baseline =  val;
    
    % Poisson rate range/values    
    case 'poisson_bb_rg'
        % Broadband poisson rate range (arbitrary units).
        NS.params.poisson_bb_rg =  val;
        
    case 'poisson_g_val'
        % Gamma poisson rate range (arbitrary units)
        NS.params.poisson_g_val =  val;
        
    case 'poisson_a_rg'
        % Alpha poisson rate range (arbitrary units)
        NS.params.poisson_a_rg =  val;
    
    % coherence ranges
    case 'gamma_coh_rg'
        % Coherence of gamma neurons within gamma band
        NS.params.gamma_coh =  val;

    % values for each stimulus
    case 'poisson_bb'
        % Poisson rates for broadband signal for each unique condition/stimulus
        % type (1 x num conditions)
        NS.params.poisson_bb = val;
    
    case 'poisson_g'
        % Poisson rates for gamma signal for each unique condition/stimulus
        % type (1 x num conditions)
        NS.params.poisson_g = val;
    
    case 'poisson_a'
        % Poisson rates for gamma signal for each unique condition/stimulus
        % type (1 x num conditions)
        NS.params.poisson_a = val;

    case 'coherence_g'
        % Poisson rates for gamma signal for each unique condition/stimulus
        % type (1 x num conditions)
        NS.params.coherence_g = val;

    % ---------------------------
    % -- trial variables --------
    % ---------------------------    
    case 'condition_num'
        % Set the condition number (equivalent to stimulus class) for each
        % trial
        NS.trial.condition_num =  val;
    case 'poisson_rate_bb'
        % Set the Poisson rate for the broadband inputs for each trial
        NS.trial.poisson_rate_bb =  val;
    case 'poisson_rate_g'
        % Set the Poisson rate for the gamma inputs for each trial
        NS.trial.poisson_rate_g =  val;
    case 'poisson_rate_a'
        % Set the Poisson rate for the alpha inputs for each trial
        NS.trial.poisson_rate_a =  val;
    case 'coherence_rate_g'
        % Set the Poisson rate for the alpha inputs for each trial
        NS.trial.coherence_rate_g =  val;

    % ---------------------------
    % -- data -------------------
    % ---------------------------
    case 'ts'
        % Set the time series data field, 3D or 4D: 
        %   time x neurons x trials x experiments
        NS.data.ts = val;
    
    case 'bb_input'
        NS.data.bb_inputs = val;
        
    case 'g_input'
        NS.data.g_inputs = val;
        
    case 'a_input'
        NS.data.a_inputs = val;
        
    case 'bb'
        % Set the ecog broadband measures (num_trials x num_experiments)
        NS.data.bb = val; 
    case 'bold'
        % Set the BOLD measures (num_trials x num_experiments)
        NS.data.bold = val; 
    case 'lfp'
        % Set the ecog lfp measures (num_trials x num_experiments)
        NS.data.lfp = val; 
    case 'gamma'
        % Set the ecog gamma measures (num_trials x num_experiments)
        NS.data.gamma = val; 
    case 'alpha'
        % Set the ecog alpha measures (num_trials x num_experiments)
        NS.data.alpha = val;         
    % ---------------------------
    % -- statistics -------------
    % ---------------------------
    case 'r2'
        % Correlation values across trial types, requires a specification
        % of which variables are correlated in varargin{1}, e.g., 
        % NS = ns_set(NS, 'r2', r2data, 'bold_bb');
        NS.stats.(varargin{1}) = val;         
     
    otherwise
        error('Unknown parameter %s', param);
end

return
function NS = ns_set(NS, param, val, varargin)
% Set parameters for neural simulation
% NS = ns_set(NS, param, val, varargin)

% Big switch
switch lower(param)        
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

    % Poisson rates for each unique condition/stimulus type (1 x num conditions)
    case 'poisson_bb'
        NS.params.poisson_bb = val;
    case 'poisson_g'
        NS.params.poisson_g = val;
    case 'poisson_a'
        NS.params.poisson_a = val;

    % Coherence for each unique condition/stimulus type (1 x num conditions)
    case 'coherence_bb'
        NS.params.coherence_bb = val;
    case 'coherence_g'
        NS.params.coherence_g = val;
    case 'coherence_a'
        NS.params.coherence_a = val;

    % ---------------------------
    % -- trial variables --------
    % ---------------------------    
    case 'condition_num'
        % Set the condition number (equivalent to stimulus class) for each
        % trial
        NS.trial.condition_num =  val;
        
    % Set the Poisson rate for each trial
    case 'trial_poisson_bb' % broadband
        NS.trial.poisson_bb =  val;
    case 'trial_poisson_g' % gamma
        NS.trial.poisson_g =  val;
    case 'trial_poisson_a' % alpha
        NS.trial.poisson_a =  val;
        
    % Set the Coherence for each trial
    case 'trial_coherence_bb' % broadband
        NS.trial.coherence_bb =  val;
    case 'trial_coherence_g' % gamma
        NS.trial.coherence_g =  val;
    case 'trial_coherence_a' % alpha
        NS.trial.coherence_a =  val;

    % ---------------------------
    % -- data -------------------
    % ---------------------------
    case 'ts'
        % Set the time series data field, 3D or 4D: 
        %   time x neurons x trials x experiments
        NS.data.ts = val;
    
    case 'lfp_spectra'
        % Set the time series data field, 3D or 4D: 
        %   time x neurons x trials x experiments
        NS.data.lfp_spectra = val;
    case 'f'
        % Set the time series data field, 3D or 4D: 
        %   time x neurons x trials x experiments
        NS.data.f = val;

    case 'bb_input'
        NS.data.bb_inputs = val;
        
    case 'g_input'
        NS.data.g_inputs = val;
        
    case 'a_input'
        NS.data.a_inputs = val;
        
    case 'bold'
        % Set the BOLD measures (num_trials)
        NS.data.bold = val; 
    case 'lfp'
        % Set the ecog lfp measures (num_trials)
        NS.data.lfp = val; 

    % Set the bootstrapped LFP values all trials (num_conditions x num_bootstraps X freq)
    case 'lfp_spectra_bs' % broadband all trials
        NS.data.lfp_spectra_bs = val; 

    % Set the bootstrapped BOLD values all trials (num_conditions x num_bootstraps)
    case 'bold_bs' % broadband all trials
        NS.data.bold_bs = val; 
    case 'bold_bs_even' % broadband all trials
        NS.data.bold_bs_even = val; 
    case 'bold_bs_odd' % broadband all trials
        NS.data.bold_bs_odd = val; 

    % Set the extracted ecog values all trials (num_conditions x num_bootstraps)
    case 'bb' % broadband all trials
        NS.data.bb = val; 
    case 'gamma' % gamma all trials
        NS.data.gamma = val; 
    case 'alpha' % alpha all trials
        NS.data.alpha = val;     

    % Set the extracted ecog values even trials (num_conditions x num_bootstraps)
    case 'bb_even' % broadband all trials
        NS.data.bb_even = val; 
    case 'gamma_even' % gamma all trials
        NS.data.gamma_even = val; 
    case 'alpha_even' % alpha all trials
        NS.data.alpha_even = val;     
        
    % Set the extracted ecog values odd trials (num_conditions x num_bootstraps)
    case 'bb_odd' % broadband all trials
        NS.data.bb_odd = val; 
    case 'gamma_odd' % gamma all trials
        NS.data.gamma_odd = val; 
    case 'alpha_odd' % alpha all trials
        NS.data.alpha_odd = val;     

    % ---------------------------
    % -- statistics -------------
    % ---------------------------
    case 'stats'
        % Correlation values across trial types, requires a specification
        % of which variables are correlated in varargin{1}, e.g., 
        % NS = ns_set(NS, 'r2', r2data, 'bold_bb');
        NS.stats = val;         
     
    otherwise
        error('Unknown parameter %s', param);
end

return
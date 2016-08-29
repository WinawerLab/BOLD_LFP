function val = ns_get(NS, param, varargin)
% Get parameters for neural simulation
% val = ns_get(NS, param, varargin)

% Big switch
switch lower(param)    
    case 'save_inputs'
        % save alpha, gamma and broadband inputs before leaky integrator
        val = NS.params.save_inputs;
    
    case 'trials_save_inputs'
        % save alpha, gamma and broadband inputs before leaky integrator
        val = NS.params.trials_save_inputs;

    case 'num_neurons'
        % The number of neurons influences how much the LFP power can deviate from
        % BOLD. This is because the sum of random signals grows with sqrt(n), and
        % synchronous signals grow with n. The larger the value of n, the more the
        % LFP signal (but not the BOLD signal!) will be dominated by the
        % synchronous signals.
        val = NS.params.num_neurons;
        
    case 'num_conditions'
        % Number of stimulus conditions. Each condition is defined by a
        % broadband and a gamma amplitude.
        val = NS.params.num_conditions;
                
    case 'num_averages'
        % Number of trials with the same condition in each experiment
        val = NS.params.num_averages;
        
    case 'num_experiments'
        % Number of experiments. (> 1 only if we want to compute a distribution of
        % R2 values)
        val = NS.params.num_experiments;
        
    case 'dt'
        % One second per trial, ms sampling
        val = NS.params.dt;
        
    case 'trial_length'
        val = NS.params.trial_length;
        
    %%%% Frequencies for band-limited responses (gamma and alpha)       
    case 'gamma_range' % Frequency range that define gamma band response (Hz)
        val = NS.params.gamma_range;
    case 'alpha_range' % Frequency range that define alpha band response (Hz)
        val = NS.params.alpha_range;
        
    % Baseline Poisson rate (arbitrary units)
    case 'poisson_baseline' 
        val = NS.params.poisson_baseline;
        
    % Poisson rates for each unique condition/stimulus type (1 x num conditions)
    case 'poisson_bb'
        val =  NS.params.poisson_bb;
    case 'poisson_g'
        val = NS.params.poisson_g;
    case 'poisson_a'
        val = NS.params.poisson_a;
    
    % Coherence rates for each unique condition/stimulus type (1 x num conditions)
    case 'coherence_bb'
        val =  NS.params.coherence_bb;
    case 'coherence_g'
        val =  NS.params.coherence_g;
    case 'coherence_a'
        val =  NS.params.coherence_a;

    % -------------------------------------------------------------------
    % Derived parameters. These can be returned by get but cannot be set
    % -------------------------------------------------------------------

    case 'num_trials'
        % Number of trials per experiment
        val = ns_get(NS, 'num_conditions') * ns_get(NS, 'num_averages');

    case 'baseline_trials'
        % Index to baseline trials (at the start of the experiment)
        val = [1:ns_get(NS,'num_averages')];

    case 't'
        % Time vector for one trial (in seconds)
        dt  = ns_get(NS, 'dt'); 
        val = dt:dt:ns_get(NS, 'trial_length');
                 
    % get frequencies for bb, gamma, alpha
    case 'freq_bb'
        f = ns_get(NS, 'f'); 
        val = f(f>90 & f<200);
    case 'freq_gamma'
        f = ns_get(NS, 'f'); 
        gamma_range = ns_get(NS, 'gamma_range');
        val = f(f > gamma_range(1) & f < gamma_range(2));
    case 'freq_alpha'
        % see 'freq_gamma'
        f = ns_get(NS, 'f'); 
        alpha_range = ns_get(NS, 'alpha_range');
        val = f(f > alpha_range(1) & f < alpha_range(2));

    case 'gamma_filter'
        % Band-pass Butterworth filter for gamma response
        dt = ns_get(NS, 'dt');  gamma_range = ns_get(NS, 'gamma_range');
        val = designfilt('bandpassiir','FilterOrder',10, ...
            'HalfPowerFrequency1',gamma_range(1),'HalfPowerFrequency2',gamma_range(2),...
            'SampleRate',1/dt,'DesignMethod','butter');

    case 'alpha_filter'
        % Band-pass Butterworth filter for alpha response
        dt = ns_get(NS, 'dt');  alpha_range = ns_get(NS, 'alpha_range');
        val = designfilt('bandpassiir','FilterOrder',10, ...
            'HalfPowerFrequency1',alpha_range(1),'HalfPowerFrequency2',alpha_range(2),...
            'SampleRate',1/dt,'DesignMethod','butter');
        
    % ---------------------------
    % -- trial variables --------
    % ---------------------------    
    case 'condition_num'
        % Get the condition number (equivalent to stimulus class) for each
        % trial
        val = NS.trial.condition_num;
        
    % Get the Poisson rate for the inputs for each trial
    case 'trial_poisson_bb'
        val = NS.trial.poisson_bb;
    case 'trial_poisson_g'
        val = NS.trial.poisson_g;
    case 'trial_poisson_a'
        val = NS.trial.poisson_a;

    % Get the coherence for the inputs for each trial
    case 'trial_coherence_bb'
        val = NS.trial.coherence_bb;
    case 'trial_coherence_g'
        val = NS.trial.coherence_g;
    case 'trial_coherence_a'
        val = NS.trial.coherence_a;
        
    % ---------------------------
    % -- data -------------------
    % ---------------------------    
    case 'ts'
        % Get the simulated time series (3D or 4D array, time x neuron x trial x experiment)
        val = NS.data.ts;
                
    case 'bold'
        % Get the bold measures (num_trials x num_experiments)
        val = NS.data.bold; 

    case 'lfp'
        % Get the ecog lfp measures (num_trials x num_experiments)
        val = NS.data.lfp; 

    case 'lfp_spectra'
        % Get the ecog lfp measures (num_trials x num_experiments)
        val = NS.data.lfp_spectra; 
    case 'f'
        val = NS.data.f; 
        
    % get the bootstrapped LFP measures (num_conditions x num_bootstraps x frequency):
    case 'lfp_spectra_bs' 
        val = NS.data.lfp_spectra_bs; 

    % get the bootstrapped BOLD measures (num_conditions x num_bootstraps):
    case 'bold_bs' 
        val = NS.data.bold_bs; 
    case 'bold_bs_even' 
        val = NS.data.bold_bs_even; 
    case 'bold_bs_odd' 
        val = NS.data.bold_bs_odd; 

    % get the ECoG measures (num_conditions x num_bootstraps):
    case 'bb' % broadband
        val = NS.data.bb; 
    case 'gamma' % gamma
        val = NS.data.gamma; 
    case 'alpha' % alpha
        val = NS.data.alpha; 
                
    % get the ECoG measures (num_conditions x num_bootstraps):
    case 'bb_even' % broadband
        val = NS.data.bb_even; 
    case 'gamma_even' % gamma
        val = NS.data.gamma_even; 
    case 'alpha_even' % alpha
        val = NS.data.alpha_even; 

    % get the ECoG measures (num_conditions x num_bootstraps):
    case 'bb_odd' % broadband
        val = NS.data.bb_odd; 
    case 'gamma_odd' % gamma
        val = NS.data.gamma_odd; 
    case 'alpha_odd' % alpha
        val = NS.data.alpha_odd; 

    otherwise
        error('Unknown parameter %s', param);
end

return
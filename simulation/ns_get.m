function val = ns_get(NS, param, varargin)
% Get parameters for neural simulation
% val = ns_get(NS, param, varargin)

% Big switch
switch lower(param)
    case 'simulate_method'
        % Two ways to simulate the neural time series:
        %   'SLOW': leaky integrator simulation
        %   'FAST': random walk
        val = NS.params.simulate_method;
    
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
                
    case 'gamma_range'
        % Band limits of spectrum that define gamma band response (Hz)
        val = NS.params.gamma_range;
        
    case 'alpha_range'
        % Band limits of spectrum that define alpha band response (Hz)
        val = NS.params.alpha_range;
        
    case 'poisson_baseline'
        % Baseline Poisson rate (arbitrary units)
        val = NS.params.poisson_baseline;
        
    % rates    
    case 'poisson_bb_rg'
        % Broadband poisson rate range (arbitrary units).
        val = NS.params.poisson_bb_rg;
                
    case 'poisson_g_val'
        % Broadband poisson rate range (arbitrary units).
        val = NS.params.poisson_g_val;
        
    case 'poisson_a_rg'
        % Broadband poisson rate range (arbitrary units).
        val = NS.params.poisson_a_rg;

    % coherences
    case 'gamma_coh_rg'
        % Coherence of gamma neurons within gamma band ([0 1])
        val =  NS.params.gamma_coh;

    % rates for each condition
    case 'poisson_bb'
        % Poisson rates for broadband signal for each unique condition/stimulus
        % type (1 x num conditions)
       val =  NS.params.poisson_bb;

    case 'poisson_g'
        % Poisson rates for gamma signal for each unique condition/stimulus
        % type (1 x num conditions)
        val = NS.params.poisson_g;
    
    case 'poisson_a'
        % Poisson rates for alpha signal for each unique condition/stimulus
        % type (1 x num conditions)
        val = NS.params.poisson_a;
    
    % coherences for each condition
    case 'coherence_g'
        % Coherence rates for gamma signal for each unique condition/stimulus
        % type (1 x num conditions)
        val =  NS.params.coherence_g;

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
        
    case 'f'
        % Temporal frequencies for one trial (in Hz)
        t = ns_get(NS, 't');
        val = (0:length(t)-1)/max(t);
                 
    case 'freq_gamma'
        % Frequency ranges for computing broadband power and gamma power are
        % complementary. Also we remove DC and frequencies above 200 Hz. For the
        % broadband calculation, we exclude frequencies that are higher than 75% of
        % the low frequency end of the gamma range, or greater than 125% of the
        % high frequency range. The 75% and 125% provide a little padding, since
        % the filter which produces the gamma signal is not hard edged, but gently
        % rolls off into nearyby frequencies. We would like to exclude these nearby
        % frequencies from the broadband calucation, because they will partially
        % reflect the gamma amplitude.
        f = ns_get(NS, 'f'); gamma_range = ns_get(NS, 'gamma_range');
        val = f(f > gamma_range(1) & f < gamma_range(2));

    case 'freq_bb'
        % see 'freq_gamma'
        f = ns_get(NS, 'f'); 
        gamma_range = ns_get(NS, 'gamma_range');
        alpha_range = ns_get(NS, 'alpha_range');
%         val = f(...
%             ((f<gamma_range(1)*.75) | (f>gamma_range(2)*1.25)) & ... % exclude gamma
%         ((f<alpha_range(1)*.5) | (f>alpha_range(2)*1.5))  ... % exclude alpha
%          & f<200 & f>6);
        val = f(f>90 & f<200);
    case 'freq_alpha'
        % see 'freq_gamma'
        f = ns_get(NS, 'f'); alpha_range = ns_get(NS, 'alpha_range');
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
        val = designfilt('bandpassiir','FilterOrder',30, ...
            'HalfPowerFrequency1',alpha_range(1),'HalfPowerFrequency2',alpha_range(2),...
            'SampleRate',1/dt,'DesignMethod','butter');
    
    case 'alpha_filter_delay'
        % Delay, in samples, caused by alpha filter
        alpha_filter = ns_get(NS, 'alpha_filter');  
        sample_rate  = 1/ns_get(NS, 'dt');
        alpha_range  = ns_get(NS, 'alpha_range');
        [delays, w] = grpdelay(alpha_filter); 
        f = w/pi * sample_rate;
        val = mean(delays(f>= alpha_range(1) & f < alpha_range(2)));
        val = round(val);
        
    case 'lowpass_filter'
        % LOW-pass Butterworth filter for alpha envelope
        dt = ns_get(NS, 'dt');  
        val = designfilt('lowpassiir','FilterOrder',10, ...
            'HalfPowerFrequency',15,...
            'SampleRate',1/dt,'DesignMethod','butter');

    % ---------------------------
    % -- trial variables --------
    % ---------------------------    
    case 'condition_num'
        % Get the condition number (equivalent to stimulus class) for each
        % trial
        val = NS.trial.condition_num;
    case 'poisson_rate_bb'
        % Get the Poisson rate for the broadband inputs for each trial
        val = NS.trial.poisson_rate_bb;
    case 'poisson_rate_g'
        % Get the Poisson rate for the broadband inputs for each trial
        val = NS.trial.poisson_rate_g;
    case 'poisson_rate_a'
        % Get the Poisson rate for the broadband inputs for each trial
        val = NS.trial.poisson_rate_a;
    case 'coherence_rate_g'
        % Get the Poisson rate for the gamma inputs for each trial
        val = NS.trial.coherence_rate_g;
        
    % ---------------------------
    % -- data -------------------
    % ---------------------------    
    case 'ts'
        % Get the simulated time series (3D or 4D array, time x neuron x trial x experiment)
        val = NS.data.ts;

    case 'bb'
        % Get the ecog broadband measures (num_trials x num_experiments)
        val = NS.data.bb; 
                
    case 'lfp'
        % Get the ecog lfp measures (num_trials x num_experiments)
        val = NS.data.lfp; 
        
    case 'gamma'
        % Get the ecog gamma measures (num_trials x num_experiments)
        val = NS.data.gamma; 
        
    case 'alpha'
        % Get the ecog alpha measures (num_trials x num_experiments)
        val = NS.data.alpha; 

    case 'bold'
        % Get the bold measures (num_trials x num_experiments)
        val = NS.data.bold; 

    case 'power'
        % Compute the population power spectra (averaging over neurons)
        % Power is frequency x trial x experiment
        ts = ns_get(NS, 'ts');
        ts_for_fft = squeeze(mean(ts,2));
        
        % fft power
        val = ns_fftpower(ts_for_fft);

    case 'fitbroadbandgamma'
        % compute broadband and gamma estimes
        data_power      = ns_get(NS, 'power'); 
        f               = ns_get(NS,'f');
        conditions      = ns_get(NS, 'condition_num');
        baseline_trials = find(conditions==0);

        f_use4fit       = 30:200;

        data_base       = mean(data_power(:,baseline_trials),2);
        out_weights     = zeros(size(data_power,2),2);
        for k=1:size(data_power,2)
            data_fit        = data_power(:,k);
            [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
                fit_gammadata(f,f_use4fit,data_base,data_fit);
            out_weights(k,:) = [w_pwr w_gauss];
        end
        val = out_weights;
        
    case 'power_law'
        % derive power law parameters for baseline
        power = ns_get(NS, 'power'); freq_bb = ns_get(NS, 'freq_bb');
        conditions = ns_get(NS, 'condition_num');
        for sim_number = 1:ns_get(NS, 'num_experiments');
            data_power_baseline = mean(log(power(freq_bb,conditions == 0, sim_number)),2);
            power_law_coefficients = polyfit(log(freq_bb)',data_power_baseline ,1);
            power_law(:,sim_number) = polyval(power_law_coefficients,log(freq_bb)');
        end
        val = power_law;
        
    otherwise
        error('Unknown parameter %s', param);
end

return
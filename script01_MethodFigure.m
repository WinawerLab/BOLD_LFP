

%% SET THE PARAMETERS

% Set default parameters.
NS = neural_sim_defaults; disp(NS.params)

% Change these to change the simulation. 
NS = ns_set(NS, 'num_neurons', 200); 
NS = ns_set(NS, 'poisson_baseline', .5); 
NS = ns_set(NS, 'poisson_bb_rg', [0 .3]); 
NS = ns_set(NS, 'poisson_g_val', .5);
NS = ns_set(NS, 'poisson_a_val', .1); 
NS = ns_set(NS, 'alpha_coh_rg', [0 .5]); 
NS = ns_set(NS, 'gamma_coh_rg', [0 .3]); 

NS = ns_set(NS,'num_conditions',2);
NS = ns_set(NS,'num_averages',1);

% Assign expected values of broadband and gamma levels for each stimulus class and each trial
NS = ns_make_trial_struct(NS); disp(NS.trial)

% overwrite settings simulation for one set
NS.params.poisson_bb_rg = [0 .3]; 
NS.params.alpha_coh = [.5 0];
NS.params.coherence_a = [.5 0];
NS.params.coherence_g = [0 .3];
NS.trial.coherence_rate_a = NS.params.coherence_a.^2;
NS.trial.coherence_rate_g = NS.params.coherence_g.^2;

%%
cond_color={[0 0 0],[0 .6 0]};

% Get variables from NS struct before simulating
t                = ns_get(NS, 't');
num_trials       = ns_get(NS, 'num_trials');
num_neurons      = ns_get(NS, 'num_neurons');
poisson_baseline = ns_get(NS, 'poisson_baseline');
poisson_rate_bb  = NS.params.poisson_bb;
poisson_rate_g   = ns_get(NS, 'poisson_rate_g');
poisson_rate_a   = ns_get(NS, 'poisson_rate_a');
coherence_rate_g = NS.params.poisson_g;
coherence_rate_a = NS.params.coherence_a;
gamma_filter     = ns_get(NS, 'gamma_filter');
alpha_filter     = ns_get(NS, 'alpha_filter');
alpha_range      = ns_get(NS, 'alpha_range');
lowpass_filter   = ns_get(NS, 'lowpass_filter');
length_zero_pad  = length(t);

%% INITIATE SIMULATION STRUCTURES

% Initialize the time series array, which will hold data for all neurons at
% all time points in all trials across all experiments
ts  = zeros(length(t), ns_get(NS, 'num_neurons'), num_trials, ...
    ns_get(NS, 'num_experiments'));

% Initialize variables. These variables will store one time
% series per tiral per neuron in the broadband pool (ts_integrated) and
% one time series per trial per neuron in the gamma pool (ts_g)
ts_integrated = zeros(length(t), num_neurons, num_trials);

% generate the simulated data: time x neuron x trial
tau     = 0.010;            % time constant of leaky integrator (seconds)
dt      = ns_get(NS, 'dt'); % step size for simulation, in seconds
t       = ns_get(NS,'t');

%% SIMULATE AND PLOT

t_idx_plot=[1:1000];
input_dc = .25;
% generate data for one trial at a time
for ii = 1:num_trials
    figure('Position',[0 0 1200 400])


    %%%%% Broadband inputs
    this_rate       = poisson_rate_bb(ii) + poisson_baseline;
    summed_inputs   = randn(length(t), num_neurons)*this_rate + input_dc;
    
    % for the first trial, plot the time-series
        subplot(4,6,1),
        plot(t(t_idx_plot),summed_inputs(t_idx_plot,1),'Color',cond_color{ii});
        title('bb inputs')
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,7),
        plot(t(t_idx_plot),summed_inputs(t_idx_plot,2),'Color',cond_color{ii});
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,13),
        plot(t(t_idx_plot),summed_inputs(t_idx_plot,3),'Color',cond_color{ii});
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])

    
    %%%%% Gamma inputs
    mu              = zeros(1,num_neurons);
    sigma           = eye(num_neurons) + (1-eye(num_neurons))* coherence_rate_g(ii);
    gamma_inputs    = mvnrnd(mu,sigma,length(t));
    gamma_inputs    = padarray(gamma_inputs, [length_zero_pad 0], 0, 'both'); % zero pad 
    gamma_inputs    = poisson_rate_g(ii)*filtfilt(gamma_filter, gamma_inputs);
    gamma_inputs    = gamma_inputs(length(t)+1:2*length(t),:); % remove zero pad 
    
        subplot(4,6,2),
        plot(t(t_idx_plot),gamma_inputs(t_idx_plot,1),'Color',cond_color{ii});
        title('gamma inputs')
        ylim([-.5 .5]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,8),
        plot(t(t_idx_plot),gamma_inputs(t_idx_plot,2),'Color',cond_color{ii});
        ylim([-.5 .5]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,14),
        plot(t(t_idx_plot),gamma_inputs(t_idx_plot,3),'Color',cond_color{ii});
        ylim([-.5 .5]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])

    %%%%% Alpha inputs with pulses
    alpha_pulses = zeros(length(t), num_neurons);
    next_pulse = 0; 
    wait_time = sort(round(1./alpha_range /dt));

   while next_pulse < length(t)
        this_pulse = next_pulse + randi(wait_time, [1 num_neurons]);
        this_pulse(this_pulse > length(t)) = NaN;
        inds = sub2ind(size(alpha_pulses), this_pulse, 1:num_neurons);
        inds = inds(isfinite(inds));
        alpha_pulses(inds) = coherence_rate_a(ii);
        next_pulse = round(mean(this_pulse));
    end
    h = exp(-(t-.075).^2/(2*.02^2));
    alpha_inputs = -conv2(alpha_pulses, h', 'full');
    [~,max_ind] = max(h);
    
    % get the peak of the response at the time of the pulse: 
    alpha_inputs = alpha_inputs(max_ind:max_ind+length(t)-1,:);
    % alpha_inputs = filtfilt(lowpass_filter, alpha_inputs); % lowpass to reduce harmonics

        subplot(4,6,3),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,1),'Color',cond_color{ii});
        title('alpha inputs')
        ylim([-.6 1]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,9),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,2),'Color',cond_color{ii});
        ylim([-.6 1]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,15),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,3),'Color',cond_color{ii});
        ylim([-.6 1]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
    
    
    % combine broadband, gamma and alpha
    summed_inputs = summed_inputs + gamma_inputs + alpha_inputs;
    

    
        subplot(4,6,4),
        plot(t(t_idx_plot),summed_inputs(t_idx_plot,1),'Color',cond_color{ii});
        title('bb+a+g')
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,10),
        plot(t(t_idx_plot),summed_inputs(t_idx_plot,2),'Color',cond_color{ii});
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,16),
        plot(t(t_idx_plot),summed_inputs(t_idx_plot,3),'Color',cond_color{ii});
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
    
    
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

    
        subplot(4,6,5),
        plot(t(t_idx_plot),ts_integrated(t_idx_plot,1,ii),'Color',cond_color{ii});
        title(['post leakyintegr'])
        axis tight 
        subplot(4,6,11),
        plot(t(t_idx_plot),ts_integrated(t_idx_plot,2,ii),'Color',cond_color{ii});
        axis tight 
        subplot(4,6,17),
        plot(t(t_idx_plot),ts_integrated(t_idx_plot,3,ii),'Color',cond_color{ii});
        axis tight 
        
        subplot(4,6,23)
        plot(t(t_idx_plot),sum(ts_integrated(t_idx_plot,:,ii),2),'Color',cond_color{ii});
        
        % plot BOLD in bar
        subplot(4,6,6)
        title('sum(x^2)')
        barh(1,mean(ts_integrated(:,1,ii).^2),.3,'FaceColor',[.5 .5 .5])
        xlim([0 1])
        subplot(4,6,12)
        barh(1,mean(ts_integrated(:,2,ii).^2),.3,'FaceColor',[.5 .5 .5])
        xlim([0 1])
        subplot(4,6,18)
        barh(1,mean(ts_integrated(:,3,ii).^2),.3,'FaceColor',[.5 .5 .5])
        xlim([0 1])
        
    for k=1:24
        subplot(4,6,k)
        if ~ismember(k,[6 12 18])
            box off
            set(gca,'XTick',[],'XColor',[1 1 1])
        else
            box off
            set(gca,'YTick',[],'YColor',[1 1 1])
        end
    end
    
%     set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['../neural_sim_output/figures/simTS2_tstable_trial' int2str(ii)])
%     print('-depsc','-r300',['../neural_sim_output/figures/simTS2_tstable_trial' int2str(ii)])
end

% subtract baseline mean
% baseline_ts = ts_integrated(:,:,ns_get(NS,'baseline_trials'));
% ts_integrated = ts_integrated - mean(baseline_ts(:));

ts = ts_integrated;
% Store the time series in the NS struct
NS = ns_set(NS, 'ts', ts);

%%
figure('Position',[0 0 400 200])
subplot(2,2,1),hold on
plot(t,sum(ts_integrated(:,:,1),2),'Color',cond_color{1});
plot(t,sum(ts_integrated(:,:,2),2),'Color',cond_color{2});
% plot(t(t_idx_plot),sum(ts_integrated(t_idx_plot,:,1),2),'Color',cond_color{1});
% plot(t(t_idx_plot),sum(ts_integrated(t_idx_plot,:,2),2),'Color',cond_color{2});
title('sum across neurons')
axis tight 

subplot(2,2,2),hold on
barh(1,sum(mean(ts_integrated(:,:,1).^2),2),.3,'FaceColor',cond_color{1})
barh(2,sum(mean(ts_integrated(:,:,2).^2),2),.3,'FaceColor',cond_color{2})
title('sum across neurons of the sum(x^2)')

subplot(2,2,3),hold on

for k=1:2
    [pxx,f]=pwelch(sum(ts(:,:,k),2),250,0,1000,1000);
    plot(f,pxx,'Color',cond_color{k})
    axis tight
    set(gca, 'XScale', 'log', 'YScale', 'log','XTick',[10 100])
    xlim([5 200])
end
    
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['../neural_sim_output/figures/simTS2_tstable_trialSum'])
% print('-depsc','-r300',['../neural_sim_output/figures/simTS2_tstable_trialSum'])



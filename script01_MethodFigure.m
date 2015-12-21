

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


not yet fixed alpha coherence not clear???


% settings simulation for one set
NS.params.poisson_bb_rg = [0 .3]; 
NS.params.alpha_coh = [.5 0];

cond_color={[0 0 0],[0 1 0]};

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
% series per tiral per neuron in the broadband pool (ts_bb) and
% one time series per trial per neuron in the gamma pool (ts_g)
ts_bb = zeros(length(t), num_neurons, num_trials);

% generate the simulated data: time x neuron x trial
tau     = 0.010;            % time constant of leaky integrator (seconds)
dt      = ns_get(NS, 'dt'); % step size for simulation, in seconds
t       = ns_get(NS,'t');

%% SIMULATE AND PLOT

t_idx_plot=[100:250];

% generate data for one trial at a time
for ii = 1:num_trials
    figure('Position',[0 0 1200 400])


    %%%%% Broadband inputs
    this_rate       = poisson_rate_bb(ii) + poisson_baseline;
    bb_inputs       = randn(length(t), num_neurons)*this_rate;
    
    % for the first trial, plot the time-series
        subplot(4,6,1),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,1),'Color',cond_color{ii});
        title('bb inputs')
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,7),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,2),'Color',cond_color{ii});
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,13),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,3),'Color',cond_color{ii});
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
        ylim([-.3 .3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,9),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,2),'Color',cond_color{ii});
        ylim([-.3 .3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,15),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,3),'Color',cond_color{ii});
        ylim([-.3 .3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
    
    
    % combine broadband, gamma and alpha
    bb_inputs = bb_inputs + gamma_inputs + alpha_inputs;
    
        subplot(4,6,4),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,1),'Color',cond_color{ii});
        title('bb+a+g')
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,10),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,2),'Color',cond_color{ii});
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
        subplot(4,6,16),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,3),'Color',cond_color{ii});
        ylim([-3 3]),xlim([t(t_idx_plot(1)) t(t_idx_plot(end))])
    
    
    % Leaky integrator loop 
    for jj = 1:length(t)-1
        % rate of change in current
        dIdt = (bb_inputs(jj,:) - ts_bb(jj,:,ii)) / tau;
        % stepwise change in current
        dI = dIdt * dt;
        % current at next time point
        ts_bb(jj+1,:,ii) = ts_bb(jj,:,ii) + dI;
    end

    
        subplot(4,6,5),
        plot(t(t_idx_plot),ts_bb(t_idx_plot,1,ii),'Color',cond_color{ii});
        title(['post leakyintegr'])
        axis tight 
        subplot(4,6,11),
        plot(t(t_idx_plot),ts_bb(t_idx_plot,2,ii),'Color',cond_color{ii});
        axis tight 
        subplot(4,6,17),
        plot(t(t_idx_plot),ts_bb(t_idx_plot,3,ii),'Color',cond_color{ii});
        axis tight 
        
        subplot(4,6,23)
        plot(t(t_idx_plot),sum(ts_bb(t_idx_plot,:,ii),2),'Color',cond_color{ii});
        
        % plot BOLD in bar
        subplot(4,6,6)
        title('sum(x^2)')
        barh(1,mean(ts_bb(:,1,ii).^2),.3,'FaceColor',[.5 .5 .5])
        xlim([0 .1])
        subplot(4,6,12)
        barh(1,mean(ts_bb(:,2,ii).^2),.3,'FaceColor',[.5 .5 .5])
        xlim([0 .1])
        subplot(4,6,18)
        barh(1,mean(ts_bb(:,3,ii).^2),.3,'FaceColor',[.5 .5 .5])
        xlim([0 .1])
        
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
    
    set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['../neural_sim_output/figures/simTS2_tstable_trial' int2str(ii)])
%     print('-depsc','-r300',['../neural_sim_output/figures/simTS2_tstable_trial' int2str(ii)])
end

% subtract baseline mean
baseline_ts = ts_bb(:,:,ns_get(NS,'baseline_trials'));
ts_bb = ts_bb - mean(baseline_ts(:));

ts = ts_bb;
% Store the time series in the NS struct
NS = ns_set(NS, 'ts', ts);

%%
figure('Position',[0 0 400 200])
subplot(2,2,1),hold on
plot(t,sum(ts_bb(:,:,1),2),'Color',cond_color{1});
plot(t,sum(ts_bb(:,:,2),2),'Color',cond_color{2});
% plot(t(t_idx_plot),sum(ts_bb(t_idx_plot,:,1),2),'Color',cond_color{1});
% plot(t(t_idx_plot),sum(ts_bb(t_idx_plot,:,2),2),'Color',cond_color{2});
title('sum across neurons')
axis tight 

subplot(2,2,2),hold on
barh(1,sum(mean(ts_bb(:,:,1).^2),2),.3,'FaceColor',cond_color{1})
barh(2,sum(mean(ts_bb(:,:,2).^2),2),.3,'FaceColor',cond_color{2})
title('sum across neurons of the sum(x^2)')
axis tight 

subplot(2,2,3),hold on

for k=1:2
    [pxx,f]=pwelch(sum(ts(:,:,k),2),250,0,1000,1000);
    plot(f,pxx,'Color',cond_color{k})
    axis tight
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlim([5 200])
end
    
set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['../neural_sim_output/figures/simTS2_tstable_trialSum'])
% print('-depsc','-r300',['../neural_sim_output/figures/simTS2_tstable_trialSum'])



%% 
%% OLD CODE
%%
%% single neuron figure
%%

figure('Position',[0 0 400 250])

% PLOT SINGLE NEURONS
subplot(4,3,[1 2 4 5 7 8]),hold on
plot(t(5:150),30+ts(5:150,1,1),'Color',[.5 0 0])
plot(t(5:150),28+ts(5:150,2,1),'Color',[.5 0 0])

plot(t(155:300),26+ts(155:300,1,2),'Color',[1 0.1 0.1])
plot(t(155:300),24+ts(155:300,2,2),'Color',[1 0.1 0.1])

plot(t(305:450),22+ts(305:450,1,3),'Color',[0 0 .5])
plot(t(305:450),20+ts(305:450,2,3),'Color',[0 0 .5])

plot(t(455:600),18+ts(455:600,1,4),'Color',[.1 .1 1])
plot(t(455:600),16+ts(455:600,2,4),'Color',[.1 .1 1])

plot(t(605:750),14+ts(605:750,301,5),'Color',[0 .5 0])
plot(t(605:750),12+ts(605:750,302,5),'Color',[0 .5 0])

plot(t(755:900),10+ts(755:900,301,6),'Color',[.1 1 .1])
plot(t(755:900),8+ts(755:900,302,6),'Color',[.1 1 .1])

xlim([.005 .9])
ylim([6 32])
set(gca,'XTick',[],'YTick',[],'YTickLabel',[])

% PLOT SUM
subplot(4,3,[10 11]),hold on
plot(t(5:150),sum(ts(5:150,1:300,1),2),'Color',[.5 0 0])
plot(t(155:300),sum(ts(155:300,1:300,2),2),'Color',[1 0.1 0.1])
plot(t(305:450),sum(ts(305:450,1:300,3),2),'Color',[0 0 .5])
plot(t(455:600),sum(ts(455:600,1:300,4),2),'Color',[.1 .1 1])
plot(t(605:750),sum(ts(605:750,301:600,5),2),'Color',[0 .5 0])
plot(t(755:900),sum(ts(755:900,301:600,6),2),'Color',[.1 1 .1])
xlim([.005 .9]),ylim([-100 100])
set(gca,'XTick',[])

% PLOT VARIANCE
subplot(4,3,[3 6 9]),hold on
barh(30,mean(ts(:,1,1).^2),'FaceColor',[.5 0 0])
barh(28,mean(ts(:,2,1).^2),'FaceColor',[.5 0 0])
barh(26,mean(ts(:,1,2).^2),'FaceColor',[1 0.1 0.1])
barh(24,mean(ts(:,2,2).^2),'FaceColor',[1 0.1 0.1])
barh(22,mean(ts(:,1,3).^2),'FaceColor',[0 0 .5])
barh(20,mean(ts(:,2,3).^2),'FaceColor',[0 0 .5])
barh(18,mean(ts(:,1,4).^2),'FaceColor',[.1 .1 1])
barh(16,mean(ts(:,2,4).^2),'FaceColor',[.1 .1 1])
barh(14,mean(ts(:,301,5).^2),'FaceColor',[0 .5 0])
barh(12,mean(ts(:,302,5).^2),'FaceColor',[0 .5 0])
barh(10,mean(ts(:,301,6).^2),'FaceColor',[.1 1 .1])
barh(8,mean(ts(:,302,6).^2),'FaceColor',[.1 1 .1])

% barh(30,var(ts(:,1,1)),'FaceColor',[.5 0 0])
% barh(28,var(ts(:,2,1)),'FaceColor',[.5 0 0])
% barh(26,var(ts(:,1,2)),'FaceColor',[1 0.1 0.1])
% barh(24,var(ts(:,2,2)),'FaceColor',[1 0.1 0.1])
% barh(22,var(ts(:,1,3)),'FaceColor',[0 0 .5])
% barh(20,var(ts(:,2,3)),'FaceColor',[0 0 .5])
% barh(18,var(ts(:,1,4)),'FaceColor',[.1 .1 1])
% barh(16,var(ts(:,2,4)),'FaceColor',[.1 .1 1])
% barh(14,var(ts(:,301,5)),'FaceColor',[0 .5 0])
% barh(12,var(ts(:,302,5)),'FaceColor',[0 .5 0])
% barh(10,var(ts(:,301,6)),'FaceColor',[.1 1 .1])
% barh(8,var(ts(:,302,6)),'FaceColor',[.1 1 .1])
ylim([6 32])
set(gca,'YTick',[8:2:30],'YTickLabel',[])

subplot(4,3,[12]),hold on
barh(6,var(sum(ts(:,1:300,1),2)),'FaceColor',[.5 0 0])
barh(5,var(sum(ts(:,1:300,2),2)),'FaceColor',[1 0.1 0.1])
barh(4,var(sum(ts(:,1:300,3),2)),'FaceColor',[0 0 .5])
barh(3,var(sum(ts(:,1:300,4),2)),'FaceColor',[.1 .1 1])
barh(2,var(sum(ts(:,301:600,5),2)),'FaceColor',[0 .5 0])
barh(1,var(sum(ts(:,301:600,6),2)),'FaceColor',[.1 1 .1])
ylim([0 7])
set(gca,'YTick',[1:6],'YTickLabel',[])

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['../neural_sim_output/figures/simulated_timeseriesSum03'])
% print('-depsc','-r300',['../neural_sim_output/figures/simulated_timeseriesSum03'])
%%
figure('Position',[0 0 400 250])
subplot(4,3,[1 2 4 5 7 8]),hold on
[pxx,f]=pwelch(sum(ts(:,1:300,1),2),250,0,1000,1000);
plot(f,pxx,'Color',[.5 0 0])
[pxx,f]=pwelch(sum(ts(:,1:300,2),2),250,0,1000,1000);
plot(f,pxx,'Color',[1 0.1 0.1])
[pxx,f]=pwelch(sum(ts(:,1:300,3),2),250,0,1000,1000);
plot(f,pxx,'Color',[0 0 .5])
[pxx,f]=pwelch(sum(ts(:,1:300,4),2),250,0,1000,1000);
plot(f,pxx,'Color',[.1 .1 1])
[pxx,f]=pwelch(sum(ts(:,301:600,5),2),250,0,1000,1000);
plot(f,pxx,'Color',[0 .5 0])
[pxx,f]=pwelch(sum(ts(:,301:600,6),2),250,0,1000,1000);
plot(f,pxx,'Color',[.1 1 .1])
axis tight
set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([5 200])
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['../neural_sim_output/figures/simulated_timeseriesPower03'])
% print('-depsc','-r300',['../neural_sim_output/figures/simulated_timeseriesPower03'])

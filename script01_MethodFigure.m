

%% SET THE PARAMETERS

% Set default parameters.
NS = neural_sim_defaults; disp(NS.params)
% Change these to change the simulation. 

NS = ns_set(NS, 'poisson_bb_rg', [0 .1]);
NS = ns_set(NS, 'poisson_g_rg', [0 1]);
NS = ns_set(NS, 'poisson_a_rg', [0 1]); 

NS = ns_set(NS,'num_conditions',3);
NS = ns_set(NS,'num_conditions',6);
NS = ns_set(NS,'num_averages',1);

% Assign expected values of broadband and gamma levels for each stimulus class and each trial
NS = ns_make_trial_struct(NS); disp(NS.trial)

NS.params.poisson_bb = [1 .1 0 0 0 0];
NS.params.poisson_a = [0 0 1 .1 0 0];
NS.params.poisson_g = [0 0 0 0 1 .1];

% Get variables from NS struct before simulating
t                = ns_get(NS, 't');
num_trials       = ns_get(NS, 'num_trials');
num_broadband    = ns_get(NS, 'num_broadband');
num_gamma        = ns_get(NS, 'num_broadband');%ns_get(NS, 'num_gamma');
num_alpha        = ns_get(NS, 'num_alpha');
poisson_rate_bb  = NS.params.poisson_bb;
poisson_baseline = ns_get(NS, 'poisson_baseline');
poisson_rate_g   = NS.params.poisson_g;
poisson_rate_a   = NS.params.poisson_a;
gamma_filter     = ns_get(NS, 'gamma_filter');
alpha_filter     = ns_get(NS, 'alpha_filter');
alpha_delay      = ns_get(NS, 'alpha_filter_delay');
envelope_filter  = ns_get(NS, 'envelope_filter');
length_zero_pad  = length(t);

%% INITIATE SIMULATION STRUCTURES

% Initialize the time series array, which will hold data for all neurons at
% all time points in all trials across all experiments
ts      = zeros(length(t), ns_get(NS, 'num_neurons'), num_trials, ...
            ns_get(NS, 'num_experiments'));

% Initialize variables. These variables will store one time
% series per tiral per neuron in the broadband pool (ts_bb) and
% one time series per trial per neuron in the gamma pool (ts_g)
ts_bb   = zeros(length(t), num_broadband, num_trials);
ts_g    = zeros(length(t), num_gamma, num_trials);

fprintf('[%s]: Simulating time series ', mfilename);
drawdots = round((1:10)/10*num_trials);

% generate the simulated data: time x neuron x trial
tau     = 0.010;            % time constant of leaky integrator (seconds)
dt      = ns_get(NS, 'dt'); % step size for simulation, in seconds
t       = ns_get(NS,'t');

%% SIMULATE AND PLOT
figure('Position',[0 0 800 400])

t_idx_plot=[1:200];

% generate data for one trial at a time
for ii = 1:num_trials
    if ismember(ii,drawdots), fprintf('.'); drawnow(); end

    %%%%% Broadband inputs
    this_rate       = poisson_rate_bb(ii) + poisson_baseline;
    bb_inputs       = randn(length(t), num_broadband)*this_rate;
    
    % for the first trial, plot the time-series
    if ii==1
        subplot(4,5,1),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,1),'r');
        title('bb inputs')
        axis tight 
        subplot(4,5,6),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,2),'r');
        axis tight 
        subplot(4,5,11),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,3),'r');
        axis tight 
    end
    
    %%%%% Gamma inputs
    mu              = zeros(1,num_gamma);
    sigma           = eye(num_gamma) + (1-eye(num_gamma))* poisson_rate_g(ii) / max(poisson_rate_g);
    gamma_inputs    = mvnrnd(mu,sigma,length(t));
    gamma_inputs    = padarray(gamma_inputs, [length_zero_pad 0], 0, 'both'); % zero pad 
    gamma_inputs    = max(poisson_rate_g)*filtfilt(gamma_filter, gamma_inputs);
    gamma_inputs    = gamma_inputs(length(t)+1:2*length(t),:); % remove zero pad 
    if ii==1
        subplot(4,5,2),
        plot(t(t_idx_plot),gamma_inputs(t_idx_plot,1),'g');
        title('gamma inputs')
        axis tight 
        subplot(4,5,7),
        plot(t(t_idx_plot),gamma_inputs(t_idx_plot,2),'g');
        axis tight 
        subplot(4,5,12),
        plot(t(t_idx_plot),gamma_inputs(t_idx_plot,3),'g');
        axis tight 
    end
    
    %%%%% Alpha inputs
    mu              = zeros(1,num_broadband); % if you add offset here it would get filtered out
    sigma           = eye(num_broadband) + (1-eye(num_broadband))* poisson_rate_a(ii) / 1;
    alpha_inputs    = mvnrnd(mu,sigma,length(t));
    alpha_inputs    = padarray(alpha_inputs, [length_zero_pad 0], 0, 'both'); % zero pad 
    alpha_inputs    = 1*filtfilt(alpha_filter, alpha_inputs);       
    alpha_envelope  = abs(hilbert(alpha_inputs));
    alpha_envelope  = alpha_envelope(length(t)+1:2*length(t),:); % remove zero pad 
    alpha_inputs    = alpha_inputs(length(t)+1:2*length(t),:); % remove zero pad 
    % add 1/alpha envelope to invert the relation between power and bold
    alpha_inputs    = alpha_inputs + (1/(1+poisson_rate_a(ii)))*alpha_envelope;
    if ii==1
        subplot(4,5,3),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,1),'b');
        title('alpha inputs')
        axis tight 
        subplot(4,5,8),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,2),'b');
        axis tight 
        subplot(4,5,13),
        plot(t(t_idx_plot),alpha_inputs(t_idx_plot,3),'b');
        axis tight 
    end
    
    % combine broadband, gamma and alpha
    bb_inputs = bb_inputs + gamma_inputs + alpha_inputs;
    if ii==1
        subplot(4,5,4),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,1),'k');
        title('bb+a+g')
        axis tight 
        subplot(4,5,9),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,2),'k');
        axis tight 
        subplot(4,5,14),
        plot(t(t_idx_plot),bb_inputs(t_idx_plot,3),'k');
        axis tight 
    end
    
    % Leaky integrator loop 
    for jj = 1:length(t)-1
        % rate of change in current
        dIdt = (bb_inputs(jj,:) - ts_bb(jj,:,ii)) / tau;
        % stepwise change in current
        dI = dIdt * dt;
        % current at next time point
        ts_bb(jj+1,:,ii) = ts_bb(jj,:,ii) + dI;
    end
    if ii==1
        subplot(4,5,5),
        plot(t(t_idx_plot),ts_bb(t_idx_plot,1,1),'k');
        title('after leaky integration')
        axis tight 
        subplot(4,5,10),
        plot(t(t_idx_plot),ts_bb(t_idx_plot,2,1),'k');
        axis tight 
        subplot(4,5,15),
        plot(t(t_idx_plot),ts_bb(t_idx_plot,3,1),'k');
        axis tight 
        
    end 
    
%     if ii==1
%         set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['../neural_sim_output/figures/simulated_timeseries01'])
%         print('-depsc','-r300',['../neural_sim_output/figures/simulated_timeseries01'])
%     end
end

ts = ts_bb;

subplot(4,5,20),
plot(t(t_idx_plot),sum(ts_bb(t_idx_plot,:,1),2),'k');
title('sum across neurons')
axis tight 

fprintf('Done\n');        

% subtract baseline mean
baseline_ts = ts(:,:,ns_get(NS,'baseline_trials'));
ts = ts - mean(baseline_ts(:));

% Store the time series in the NS struct
NS = ns_set(NS, 'ts', ts);

% save(['../neural_sim_output/data/examplesetB'],'NS') % set B is nice

subplot(4,5,19),hold on

cond_color={[.5 0 0],[1 0.1 0.1],[0 0 .5],[.1 .1 1],[0 .5 0],[.1 1 .1]};

for k=1:6
    [pxx,f]=pwelch(sum(ts(:,:,k),2),250,0,1000,1000);
    plot(f,pxx,'Color',cond_color{k})
    axis tight
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlim([5 200])
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['../neural_sim_output/figures/simTS_tstable'])
print('-depsc','-r300',['../neural_sim_output/figures/simTS_tstable'])


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

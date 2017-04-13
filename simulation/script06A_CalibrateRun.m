%%
% This script generates pannels for Fig 6 from Hermes et al:
%
% Purpose: Simulate neural data - time varying membrane potentials - with a
% structured set of inputs, such that the relation between the input level
% (amplitude and coherence) of C1, C2 and C3 and the output broadband,
% gamma and alpha power can be described
%
% DH 2016


%% Simulate data for calibration:

sim_nr = 2; 
% sim_nr = 2 has input_dc = .25 and alpha_envelop*1, poisson_baseline = .3
% and alpha_synchrony  = .75

% make conditions of input -  vary gamma, alpha and broadband
% systematically: change amplitude and coherence of each
nr_conds = 10;
out = script_make_calibration_conds(nr_conds);

tic
clear ns_params
for k=1:size(out.poisson_bb,2) % number of settings

    % Set default parameters
    NS = neural_sim_defaults; %disp(NS.params)

    NS = ns_set(NS, 'save_inputs', 1);
    NS = ns_set(NS, 'num_conditions', size(out.poisson_bb,1));
    NS = ns_set(NS, 'num_neurons', 200); 
    NS = ns_set(NS, 'poisson_baseline', .3); 
    NS = ns_set(NS, 'poisson_bb',out.poisson_bb(:,k));
    NS = ns_set(NS, 'poisson_g',out.poisson_g(:,k));
    NS = ns_set(NS, 'poisson_a',out.poisson_a(:,k));
    NS = ns_set(NS, 'coherence_bb',out.coherence_bb(:,k));
    NS = ns_set(NS, 'coherence_g',out.coherence_g(:,k));
    NS = ns_set(NS, 'coherence_a',out.coherence_a(:,k));
    
    % set bb, gamma and alpha values for all trials
    NS = ns_make_trial_struct(NS);
    
    % if save_inputs choose trials to save, data can get big if saving a lot
    NS = ns_set(NS, 'trials_save_inputs',[1 length(NS.trial.poisson_bb)]);

    % Simulate. This will produce a time series for each neuron in each trial
    NS = ns_simulate_data(NS); 

    % Convert the neural time series into instrument (LFP/BOLD) measures
    NS = ns_neural2instruments(NS); %disp(NS.data)

    % Get spectra and broadband, gamma, alpha and mean from the LFP
    NS = ns_analyse_lfp(NS); %disp(NS.data)
    
    % Compute the correlations between different instrument measures 
    NS = ns_summary_statistics(NS); %disp(NS.stats)
        
    NS.data.ts = single(NS.data.ts); % to reduce size
    
    save(fullfile(BOLD_LFPRootPath, 'data', sprintf('NS_simnr%d_set%d', sim_nr,k)),'NS');
    ns_params{k} = NS.params;
    disp(['done simulation ' int2str(k) ' of ' int2str(size(out.poisson_bb,2))])
    toc
end

%% plot lookup values and create lookup vectors for gamma and alpha

clear all

sim_nr = 2; % simulation number (general settings change)

lookup = struct([]); 

param_set_vary = {'poisson_g','coherence_g','poisson_a','coherence_a','poisson_bb','coherence_bb'};

f = figure('Position',[0 0 700 250]);hold on
for param_set = 1:4
    %load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(param_set) ],'NS')
    load(fullfile(BOLD_LFPRootPath, 'data', sprintf('NS_simnr%d_set%d', sim_nr, param_set)), 'NS');

    bold_avg        = median(NS.data.bold_bs,2);
    b_avg           = median(NS.data.bb,2);
    g_avg           = median(NS.data.gamma,2);
    a_avg           = median(NS.data.alpha,2);
    a_avg           = a_avg+1.5;
    % bold_ci         = [quantile(NS.data.bold_bs,.16,2) quantile(NS.data.bold_bs,.84,2)];
    % bb_ci           = [quantile(NS.data.bb,.16,2) quantile(NS.data.bb,.84,2)];
    % gamma_ci        = [quantile(NS.data.gamma,.16,2) quantile(NS.data.gamma,.84,2)];
    % alpha_ci        = [quantile(NS.data.alpha,.16,2) quantile(NS.data.alpha,.84,2)];
    
    % get inputs
    eval([...
        'x = NS.params.' param_set_vary{param_set} ';'...
        'y = ' param_set_vary{param_set}(end) '_avg;']);

    % estimate 
    X0 = [1 1 1 0 1];
    LB = [0 0 0 0 1];
    UB = [Inf Inf Inf Inf Inf];
    my_options=optimset('Display','off','Algorithm','trust-region-reflective'); 
    [b] = lsqnonlin(@(b) ns_in2out(b,y,x,b_avg),X0,LB,UB,my_options);
    
    subplot(1,4,param_set),hold on
    % plot data
    plot(x,y,'k.','MarkerSize',20)
    % plot estimate with inputs
    plot(x,b(1) .* 10.^(-b_avg./b(5)) .* log10((b(2)+x)./b(2) + b_avg*b(3)) + b(4),'rs');
    plot(x(1:10),b(1) .* 10.^(-b_avg(1:10)./b(5)) .* log10((b(2)+x(1:10))./b(2) + b_avg(1:10)*b(3)) + b(4),'r');
    plot(x(11:20),b(1) .* 10.^(-b_avg(11:20)./b(5)) .* log10((b(2)+x(11:20))./b(2) + b_avg(11:20)*b(3)) + b(4),'r');
    plot(x(21:30),b(1) .* 10.^(-b_avg(21:30)./b(5)) .* log10((b(2)+x(21:30))./b(2) + b_avg(21:30)*b(3)) + b(4),'r');

    xlabel([param_set_vary{param_set}]),ylabel(['LFP ' param_set_vary{param_set}(end)])
    lookup = setfield(lookup,{param_set},'b',{1:length(b)},b);
    
end

set(gcf,'PaperPositionMode','auto')

print('-depsc','-r300',fullfile(BOLD_LFPRootPath, 'figures', sprintf('sim%d_calibrate_lfpvals_g_a_new', sim_nr)))


%% broadband lookup table and save all

f = figure('Position',[0 0 700 250]);hold on
for param_set = 5:6
    %load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(param_set) ],'NS')
    load(fullfile(BOLD_LFPRootPath, 'data', sprintf('NS_simnr%d_set%d', sim_nr, param_set)), 'NS');

    bold_avg    = median(NS.data.bold_bs,2);
    b_avg       = median(NS.data.bb,2);
    
    % get inputs
    eval([...
        'x = NS.params.' param_set_vary{param_set} ';'...
        'y = ' param_set_vary{param_set}(end) '_avg;']);
    
    % estimate 
    X0 = [1 1];
    LB = [0 0];
    UB = [Inf Inf];

    my_options=optimset('Display','off','Algorithm','trust-region-reflective'); 
    [b] = lsqnonlin(@(b) ns_in2out_singleparam(b,y,x),X0,LB,UB,my_options);

    subplot(1,4,param_set-4),hold on
    % plot data
    plot(x,y,'k.','MarkerSize',20)
    % plot estimate with inputs
    plot(x,b(1) .* (log10((b(2)+x)./b(2))),'rd')
    plot(x,b(1) .* (log10((b(2)+x)./b(2))),'r')
    xlabel([param_set_vary{param_set}]),ylabel(['LFP ' param_set_vary{param_set}(end)])
    
    lookup = setfield(lookup,{param_set},'b',{1:length(b)},b);
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/calibrate_lfpvals_bb_new'])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/calibrate_lfpvals_bb_new'])

%save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_lookup_table' ],'lookup')
save(fullfile(BOLD_LFPRootPath, 'data', sprintf('NS_simnr%d_lookup_table', sim_nr)), 'lookup');
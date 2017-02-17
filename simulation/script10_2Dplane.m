%% ECOG BOLD simulation two main simulations

% Purpose: Simulate neural data - time varying membrane potentials - and
% then ask whether the simulated BOLD signal and various metrics of the
% simulated field potentials are correlated.

% sim_nr = 1; 
% sim_nr = 1 has input_dc = .2 and alpha_envelop*.5, poisson_baseline = .2
% previous fitting function without offset and b1/(b2+bb)

% sim_nr = 2; 
% sim_nr = 2 has input_dc = .25 and alpha_envelop*1, poisson_baseline = .3
% and alpha_synchrony  = .75

% sim_nr = 3; 
% sim_nr = 3 has input_dc = .25 and alpha_envelop*1, poisson_baseline = .4
% and alpha_synchrony  = .75

sim_nr = 4; 

% make conditions of input:
nr_conds = 5;
out = script_make_2Dplane_conds(nr_conds);

tic
clear ns_params
for k=1:size(out.poisson_bb,2) % number of settings

    % Set default parameters
    NS = neural_sim_defaults; %disp(NS.params)

    NS = ns_set(NS, 'save_inputs', 1);
    NS = ns_set(NS, 'num_conditions', size(out.poisson_bb,1));
    NS = ns_set(NS, 'num_neurons', 200); 
    NS = ns_set(NS, 'poisson_baseline', .4); 
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

%     save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(k) '_simonly'],'NS')
    
    % Convert the neural time series into instrument (LFP/BOLD) measures
    NS = ns_neural2instruments(NS); %disp(NS.data)

    % Get spectra and broadband, gamma, alpha and mean from the LFP
    NS = ns_analyse_lfp(NS); %disp(NS.data)
    
    % Compute the correlations between different instrument measures 
    NS = ns_summary_statistics(NS); %disp(NS.stats)
        
    NS.data.ts = single(NS.data.ts); % to reduce size
    
%     save(['../data/NS_simnr' int2str(sim_nr) '_set' int2str(k) ],'NS')
    save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(k) ],'NS')
    ns_params{k} = NS.params;
    disp(['done simulation ' int2str(k) ' of ' int2str(size(out.poisson_bb,2))])
    toc
end

%% plot lookup values and create lookup vectors for gamma and alpha

clear all

sim_nr = 4; 

% get conditions of input:
nr_conds = 5; 
out = script_make_2Dplane_conds(nr_conds);

bb_in_poisson = zeros(size(out.poisson_bb,1),size(out.poisson_bb,2));
bb_in_coherence = zeros(size(out.poisson_bb,1),size(out.poisson_bb,2));
bold_avg = zeros(size(out.poisson_bb,1),size(out.poisson_bb,2));
b_avg = zeros(size(out.poisson_bb,1),size(out.poisson_bb,2));
g_avg = zeros(size(out.poisson_bb,1),size(out.poisson_bb,2));
a_avg = zeros(size(out.poisson_bb,1),size(out.poisson_bb,2));

for param_set = 1:size(out.poisson_bb,2) % number of settings
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(param_set) ],'NS')

    bb_in_poisson(:,param_set)      = NS.params.poisson_bb;
    bb_in_coherence(:,param_set)    = NS.params.coherence_bb;
    bold_avg(:,param_set)           = median(NS.data.bold_bs,2);
    b_avg(:,param_set)              = median(NS.data.bb,2);
    g_avg(:,param_set)              = median(NS.data.gamma,2);
    a_avg(:,param_set)              = median(NS.data.alpha,2);
    a_avg(:,param_set)              = a_avg(:,param_set)+1.5;
    % bold_ci         = [quantile(NS.data.bold_bs,.16,2) quantile(NS.data.bold_bs,.84,2)];
    % bb_ci           = [quantile(NS.data.bb,.16,2) quantile(NS.data.bb,.84,2)];
    % gamma_ci        = [quantile(NS.data.gamma,.16,2) quantile(NS.data.gamma,.84,2)];
    % alpha_ci        = [quantile(NS.data.alpha,.16,2) quantile(NS.data.alpha,.84,2)];
end
%%

% f = figure('Position',[0 0 700 250]);hold on
figure('Position',[0 0 700 250]);hold on
subplot(1,2,1)
plot3(bb_in_poisson,bb_in_coherence,b_avg,'k.','MarkerSize',20);
xlabel('poisson')
ylabel('coherence')
zlabel('LFP bb')

subplot(1,2,2)
plot3(bb_in_poisson,bb_in_coherence,bold_avg,'k.','MarkerSize',20);
xlabel('poisson')
ylabel('coherence')
zlabel('BOLD')

print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/2Dplane'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/2Dplane'])








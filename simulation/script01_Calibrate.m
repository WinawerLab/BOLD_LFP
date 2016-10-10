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

% make conditions of input:
nr_conds = 10;
out.poss = script_make_calibration_conds(nr_conds);

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

sim_nr = 3; % simulation number (general settings change)

lookup = struct([]); 

param_set_vary = {'poisson_g','coherence_g','poisson_a','coherence_a','poisson_bb','coherence_bb'};

f = figure('Position',[0 0 700 250]);hold on
for param_set = 1:4
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(param_set) ],'NS')

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
    X0 = [1 1 0];
    LB = [-Inf -Inf 0];
    UB = [Inf Inf Inf];
    my_options=optimset('Display','off','Algorithm','trust-region-reflective'); 
    [b] = lsqnonlin(@(b) ns_in2out(b,y,x,b_avg),X0,LB,UB,my_options);
    
    subplot(1,4,param_set),hold on
    % plot data
    plot(x,y,'k.','MarkerSize',20)
    % plot estimate with inputs
%     plot(x,b(4)./(1+b_avg) .* (b(1).*x.^2 + b(2).*x + b(3)),'ms')
    plot(x,b(1)./(1+b_avg) .* x.^b(2) + b_avg*b(3),'cs')
    xlabel([param_set_vary{param_set}]),ylabel(['LFP ' param_set_vary{param_set}(end)])
    lookup = setfield(lookup,{param_set},'b',{1:length(b)},b);
    
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/calibrate_lfpvals_g_a'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/calibrate_lfpvals_g_a'])

f = figure('Position',[0 0 700 250]);hold on
for param_set = 5:6
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set' int2str(param_set) ],'NS')

    bold_avg    = median(NS.data.bold_bs,2);
    b_avg       = median(NS.data.bb,2);
    
    % get inputs
    eval([...
        'x = NS.params.' param_set_vary{param_set} ';'...
        'y = ' param_set_vary{param_set}(end) '_avg;']);
    
    % estimate 
    X0 = [1 0];
    LB = [-Inf -Inf];
    UB = [Inf Inf];
    my_options=optimset('Display','off','Algorithm','trust-region-reflective'); 
    [b] = lsqnonlin(@(b) ns_in2out_singleparam(b,y,x),X0,LB,UB,my_options);

    subplot(1,4,param_set-4),hold on
    % plot data
    plot(x,y,'k.','MarkerSize',20)
    % plot estimate with inputs
%     plot(x,b(4)./(1+bb_avg) .* (b(1).*x.^2 + b(2).*x + b(3)),'ms')
    plot(x,b(1) .* x.^b(2),'cs')
    xlabel([param_set_vary{param_set}]),ylabel(['LFP ' param_set_vary{param_set}(end)])
    
    lookup = setfield(lookup,{param_set},'b',{1:length(b)},b);
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/calibrate_lfpvals_bb'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/calibrate_lfpvals_bb'])

save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_lookup_table' ],'lookup')


%% 

clear all
sim_nr = 2;
% load the lookup table
load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_lookup_table' ],'lookup');

% load the data
load('/Volumes/DoraBigDrive/data/visual/m-files/bold_datalikesimulation/data/boldecog_structure_final.mat');
    
lookup_combs=[...
    1 3 5
    1 4 5
    2 3 5
    2 4 5
    1 3 6
    1 4 6
    2 3 6
    2 4 6];

%%
out = script_make_calibration_conds(10);

y = [1 1 -.3]; % 1 electrode, 1 condition

bb_amp = y(1);
if bb_amp<0 % can not bb values < 0
    bb_amp = 0;
end
gamma_amp = y(2);
alpha_amp = y(3)+1.2;

k = 3; % first lookup combination
% broadband
b = lookup(lookup_combs(k,3)).b; % lookup_combs(k,3) is 5 or 6: vary level/coherence
if lookup_combs(k,3)==5 % lookup bb level, coherence fixed
    poisson_bb = nthroot(bb_amp./b(1),b(2));
    coherence_bb = out.coherence_bb(1,5);
elseif lookup_combs(k,3)==6 % lookup bb coherence, level fixed
    poisson_bb = out.poisson_bb(1,6);
    coherence_bb = nthroot(bb_amp./b(1),b(2));
end

% gamma
b = lookup(lookup_combs(k,1)).b; % lookup_combs(k,1) is 1 or 2: vary level/coherence
if lookup_combs(k,1)==1 % lookup gamma level, coherence fixed
    poisson_g = nthroot(gamma_amp*(b(2)+bb_amp)./b(1),b(3));
    coherence_g = out.coherence_g(1,1);
elseif lookup_combs(k,1)==2 % lookup gamma coherence, level fixed
    poisson_g = out.poisson_g(1,2);
    coherence_g = nthroot(gamma_amp*(b(2)+bb_amp)./b(1),b(3));
end

% alpha
b = lookup(lookup_combs(k,2)).b; % lookup_combs(k,2) is 3 or 4: vary level/coherence
if lookup_combs(k,2)==3 % lookup alpha level, coherence fixed
    poisson_a = nthroot(alpha_amp*(b(2)+bb_amp)./b(1),b(3));
    coherence_a = out.coherence_a(1,3);
elseif lookup_combs(k,2)==4 % lookup alpha coherence, level fixed
    poisson_a = out.poisson_a(1,4);
    coherence_a = nthroot(alpha_amp*(b(2)+bb_amp)./b(1),b(3));
end








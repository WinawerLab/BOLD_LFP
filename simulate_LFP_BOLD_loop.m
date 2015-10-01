
%% ECOG BOLD simulation

% Purpose: Simulate neural data - time varying membrane potentials - and
% then ask whether the simulated BOLD signal and various metrics of the
% simulated field potentials are correlated.

%%%%% vary broadband 
% sim_nr = 1;
% poisson_bb_rg_in = [0 .5; 0 .4;0 .3;0 .2;0 .1;0 0];
sim_nr = 2;
poisson_bb_rg_in = [0 .2;0 .1;0 0.08;0 0.06;0 0.04;0 0.02; 0 0];
tic
for k=1:size(poisson_bb_rg_in,1)

    % Set default parameters.
    NS = neural_sim_defaults; %disp(NS.params)
    % Change these to change the simulation. 

    NS = ns_set(NS, 'poisson_bb_rg', poisson_bb_rg_in(k,:));
    NS = ns_set(NS, 'poisson_g_rg', [0 .5]);
    NS = ns_set(NS, 'poisson_a_rg', [0 .5]); 

    % Assign expected values of broadband and gamma levels for each stimulus class and each trial
    NS = ns_make_trial_struct(NS); %disp(NS.trial)

    % Simulate. This will produce a time series for each neuron in each trial
    NS = ns_simulate_data(NS); 

    % Convert the neural time series into instrument measures
    NS = ns_neural2instruments(NS); %disp(NS.data)

    % Compute the correlations between different instrument measures 
    NS = ns_summary_statistics(NS); %disp(NS.stats)
    
    save(['../neural_sim_output/data/NS_varyBB0-0p' int2str(10*poisson_bb_rg_in(k,2)) '_nr' int2str(sim_nr)],'NS')
    disp(['done simulation ' int2str(k) ' of ' int2str(size(poisson_bb_rg_in,1))])
    toc
end

%%
clear all
% sim_nr = 1;
% poisson_bb_rg_in = [0 .5; 0 .4;0 .3;0 .2;0 .1;0 0];
sim_nr = 2;
poisson_bb_rg_in = [0 .2;0 .1;0 0.08;0 0.06;0 0.04;0 0.02; 0 0];
bold_corr=zeros(size(poisson_bb_rg_in,1),4);

for k=1:size(poisson_bb_rg_in,1)
    load(['../neural_sim_output/data/NS_varyBB0-0p' int2str(10*poisson_bb_rg_in(k,2)) '_nr' int2str(sim_nr)],'NS')
    
    alpha_avg = ns_mean_by_stimulus(NS, ns_get(NS, 'alpha'));
    bb_avg  = ns_mean_by_stimulus(NS, ns_get(NS, 'bb'));
    gamma_avg  = ns_mean_by_stimulus(NS, ns_get(NS, 'gamma'));
    bold_avg  = ns_mean_by_stimulus(NS, ns_get(NS, 'bold'));
    all_mean_data = squeeze(mean(sum(NS.data.ts(:,:,:),2),1));
    mean_avg = zeros(max(NS.trial.condition_num),1);
    
    for m=1:max(NS.trial.condition_num)+1
        mean_avg(m) = mean(all_mean_data(NS.trial.condition_num==m-1));
    end
    
    bold_corr(k,1)=corr(bb_avg,bold_avg);
    bold_corr(k,2)=corr(gamma_avg,bold_avg);
    bold_corr(k,3)=corr(alpha_avg,bold_avg);
    bold_corr(k,4)=corr(mean_avg,bold_avg);

end
%%

figure
% bar(bold_corr(:,1))
bar(bold_corr)


%% ECOG BOLD simulation two main simulations

% Purpose: Simulate neural data - time varying membrane potentials - and
% then ask whether the simulated BOLD signal and various metrics of the
% simulated field potentials are correlated.

out.poisson_bb =   [.0  .4  ;.0  .0; .0  .0]';
out.poisson_g =    [.0  .0  ;.2  .2; .0  .0]';
out.poisson_a =    [.0  .0  ;.0  .0; .1  .5]';
%%%%% vary coherence ranges across sets:
out.coherence_bb = [.0  .0  ;.0  .0; .0  .0]';
out.coherence_g =  [.0  .0  ;.0  .6; .0  .0]';
out.coherence_a =  [.0  .0  ;.0  .0; .75 .75]
% figure('Position',[0 0 800 200])

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
            
    NS.data.ts = single(NS.data.ts); % to reduce size
    
    ns_params{k} = NS.params;
    disp(['done simulation ' int2str(k) ' of ' int2str(size(out.poisson_bb,2))])
    toc

    y0 = mean(NS.data.lfp_spectra(:,NS.trial.condition_num==0),2);
    y1 = mean(NS.data.lfp_spectra(:,NS.trial.condition_num==1),2);

    subplot(1,3,k),hold on
    plot(NS.data.f,y0,'Color',[0 1 1])
    plot(NS.data.f,y1,'Color',[0 0 1])
%     plot(NS.data.f,y0,'Color',[.5 .5 .5])
%     plot(NS.data.f,y1,'Color',[0 0 0])
    set(gca, 'XScale', 'log', 'YScale', 'log','XTick',[10 100],'YTick',[10.^-2 10.^0 10.^2])
    xlabel ('Frequency'), ylabel('Power')
    xlim([1 max(NS.data.f)]); %ylim([10.^-3 10.^2]);

end
%%
for k = 1:3
    subplot(1,3,k)
%     ylim([10.^-4 10.^2+10])
    ylim([10.^-4 10.^4+10])
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/Paper_Figure3_2Windows'])
print('-dpng','-r300',['../figures/Paper_Figure3_2Windows'])


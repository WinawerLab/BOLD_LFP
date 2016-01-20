
%% ECOG BOLD simulation two main simulations

% Purpose: Simulate neural data - time varying membrane potentials - and
% then ask whether the simulated BOLD signal and various metrics of the
% simulated field potentials are correlated.

%%%%% vary broadband 
sim_nr = 1;
poisson_bb_rg_in = [0 .5;0 .2];
tic

for k=1:size(poisson_bb_rg_in,1)

    % Set default parameters.
    NS = neural_sim_defaults; %disp(NS.params)
    
    % Change these to change the simulation. 
    NS = ns_set(NS, 'save_inputs', 1);
    NS = ns_set(NS, 'num_neurons', 200); 
    NS = ns_set(NS, 'poisson_baseline', .5); 
    NS = ns_set(NS, 'poisson_bb_rg', poisson_bb_rg_in(k,:));
    NS = ns_set(NS, 'poisson_g_val', .5);
    NS = ns_set(NS, 'poisson_a_rg', [0 .3]); 
    NS = ns_set(NS, 'gamma_coh_rg', [0 1]); 
    NS = ns_set(NS, 'num_conditions',8);

    % Assign expected values of broadband and gamma levels for each stimulus class and each trial
    NS = ns_make_trial_struct(NS); 

    % if save_inputs choose trials to save, data can get big if saving a lot
    NS = ns_set(NS, 'trials_save_inputs',[1 length(NS.trial.poisson_rate_bb)]);

    % Simulate. This will produce a time series for each neuron in each trial
    NS = ns_simulate_data(NS); 

    % Convert the neural time series into instrument measures
    NS = ns_neural2instruments(NS); %disp(NS.data)

    % Compute the correlations between different instrument measures 
    NS = ns_summary_statistics(NS); %disp(NS.stats)
    
    NS.data.ts = single(NS.data.ts); % to reduce size
    
    save(['./data/NS_simnr' int2str(sim_nr) '_set' int2str(k) ],'NS')
%     save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(k) ],'NS')
    disp(['done simulation ' int2str(k) ' of ' int2str(size(poisson_bb_rg_in,1))])
    toc
end


%% get all the correlation values 
clear all
sim_nr = 1;
poisson_bb_rg_in = [0 .5;0 1];

bold_corr=zeros(size(poisson_bb_rg_in,1),7);
bold_beta=zeros(size(poisson_bb_rg_in,1),12);
beta_ind = [1 2 3 3 4 5 5 6 6 7 7 7];% which beta belongs to which model

r_all = zeros(size(poisson_bb_rg_in,1),501);
out.bold_vals = zeros(size(poisson_bb_rg_in,1),8);
out.bb_vals = zeros(size(poisson_bb_rg_in,1),8);
out.gamma_vals = zeros(size(poisson_bb_rg_in,1),8);
out.alpha_vals = zeros(size(poisson_bb_rg_in,1),8);
out.mean_vals = zeros(size(poisson_bb_rg_in,1),8);

for k=1:size(poisson_bb_rg_in,1)
    load(['./data/NS_simnr' int2str(sim_nr) '_set' int2str(k) ],'NS')
%     load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(k) ],'NS')
    
    %- to correlate bold and neurophys
    alpha_avg = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'alpha')));
    bb_avg  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bb')));
    gamma_avg  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'gamma')));
    bold_avg  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bold')));
    all_mean_data = squeeze(mean(sum(NS.data.ts(:,:,:),2),1));
    mean_avg = zeros(max(NS.trial.condition_num),1);
    for m=1:max(NS.trial.condition_num)+1
        mean_avg(m) = mean(all_mean_data(NS.trial.condition_num==m-1));
    end
    mean_avg = zscore(mean_avg);
    
    % put in output
    out.bold_vals(k,:) = bold_avg;
    out.bb_vals(k,:) = bb_avg;
    out.gamma_vals(k,:) = gamma_avg;
    out.alpha_vals(k,:) = alpha_avg;
    out.mean_vals(k,:) = mean_avg;
    
    % do stats for different models
    stats = regstats(bold_avg,bb_avg);
    bold_corr(k,1) = stats.rsquare;
    bold_beta(k,1) = stats.beta(2);

    stats = regstats(bold_avg,gamma_avg);
    bold_corr(k,2) = stats.rsquare;
    bold_beta(k,2) = stats.beta(2);
    
    stats = regstats(bold_avg,[bb_avg gamma_avg]);
    bold_corr(k,3) = stats.rsquare;
    bold_beta(k,3:4) = stats.beta(2:3);

    stats = regstats(bold_avg,alpha_avg);
    bold_corr(k,4) = stats.rsquare;
    bold_beta(k,5) = stats.beta(2);
    
    stats = regstats(bold_avg,[bb_avg alpha_avg]);
    bold_corr(k,5) = stats.rsquare;
    bold_beta(k,[6:7]) = stats.beta(2:3);
    
    stats = regstats(bold_avg,[gamma_avg alpha_avg]);
    bold_corr(k,6) = stats.rsquare;
    bold_beta(k,[8:9]) = stats.beta(2:3);

    stats = regstats(bold_avg,[bb_avg gamma_avg alpha_avg]);
    bold_corr(k,7) = stats.rsquare;
    bold_beta(k,[10:12]) = stats.beta(2:4);

    %- correlate across all frequencies
    bold_avg  = ns_mean_by_stimulus(NS, ns_get(NS, 'bold'));
    ts_for_fft = squeeze(mean(ns_get(NS, 'ts'),2));
    
    % fft settings just like in the other correlation
    fft_w  = hanning(250); % window width
    fft_ov = .5*length(fft_w); % overlap
    srate = 1/ns_get(NS,'dt');
    % initialize 
    [~,f] = pwelch(ts_for_fft(:,1),fft_w,fft_ov,srate,srate);
    pxx_all = zeros(size(ts_for_fft,2),length(f));
    % loop over trials
    for m = 1:size(ts_for_fft,2)
        pxx_all(m,:) = pwelch(ts_for_fft(:,m),fft_w,fft_ov,srate,srate);
    end
    % mean power by stimulus
    pxx_avg = zeros(NS.params.num_conditions,length(f));
    for m=1:NS.params.num_conditions
        pxx_avg(m,:) = mean(pxx_all(NS.trial.condition_num==m-1,:));
    end
    % power change
    pxx_change = bsxfun(@minus, log10(pxx_avg), log10(pxx_avg(1,:)));
    r = zeros(size(pxx_change,2),1);
    for m=1:size(pxx_change,2)
        r(m) = corr(bold_avg,pxx_change(:,m));
    end
    r_all(k,:) = r;
    
end

%%
%% figures for checks:
%%

%% Plot beta and r2 for each model separately
bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};
figure('Position',[0 0 580 200])

for s = 1:size(poisson_bb_rg_in,1)
    subplot(1,2,s),hold on
    for k=1:7
        bar(k,bold_corr(s,k),.9,'FaceColor',bar_colors{k})
    end
    xlim([0 8]),ylim([0 1])
    set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb g','a','bb a','g a','bb g a'})
    title(['Model ' int2str(s) ' R^2'])
end


figure('Position',[0 0 220 200])
labels_beta={{'bb','',''},{'','g',''},{'bb','g',''},{'','','a'},...
    {'bb','','a'},{'','g','a'},{'bb','g','a'}};
labels_index={[1],[2],[1 2],[3],[1 3],[2 3],[1 2 3]};
bb_g_a_color={[.9 .9 .9],[.6 .6 .6],[.3 .3 .3]};

for s = 1:size(poisson_bb_rg_in,1)
    for k=1:7
        xl_ind=labels_index{k};
        subplot(size(poisson_bb_rg_in,1),7,(s-1)*7+k),hold on 

        temp_beta=bold_beta(s,beta_ind==k);      
        for m=1:length(xl_ind)
            bar(xl_ind(m),temp_beta(m),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
        end

        xlim([.5 3.5]),ylim([-1 1.5])
        set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-1:.2:1.5],'YTickLabel',[])
    end
end

%% correlation across all frequencies

figure('Position',[0 0 400 120])
for s = 1:size(poisson_bb_rg_in,1)
    subplot(1,2,s),hold on
    plot(f,zeros(size(f)),'Color',[.5 .5 .5])
    plot(f,r_all(s,:),'k','LineWidth',2)
    xlabel('Frequency (Hz)'), ylabel('correlation (r)')
    xlim([0 200])
    ylim([-1 1])
end

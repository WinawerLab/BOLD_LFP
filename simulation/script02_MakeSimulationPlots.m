
%% get all the correlation values 
clear all
sim_nr      = 1;
ns_files    = dir(['./data/NS_simnr' int2str(sim_nr) '_set*']);

bold_corr=zeros(length(ns_files),7);
bold_beta=zeros(length(ns_files),12);
beta_ind = [1 2 3 3 4 5 5 6 6 7 7 7];% which beta belongs to which model

r_all = zeros(length(ns_files),501);
out.bold_vals = zeros(length(ns_files),8);
out.bb_vals = zeros(length(ns_files),8);
out.gamma_vals = zeros(length(ns_files),8);
out.alpha_vals = zeros(length(ns_files),8);
out.mean_vals = zeros(length(ns_files),8);
NSall=[];

for k=1:length(ns_files)
    load(['./data/' ns_files(k).name],'NS')
    
%     load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_nr' int2str(sim_nr) '_' int2str(k) ],'NS')
    NSall{k} = NS;
    
    % to correlate bold and neurophys
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
    % even values
    alpha_even = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'alpha'),'even'));
    bb_even  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bb'),'even'));
    gamma_even  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'gamma'),'even'));
    bold_even  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bold'),'even'));
    all_mean_data = squeeze(mean(sum(NS.data.ts(:,:,:),2),1));
    mean_even = zeros(max(NS.trial.condition_num),1);
    for m=1:max(NS.trial.condition_num)+1
        mean_even(m) = mean(all_mean_data(NS.trial.condition_num==m-1));
    end
    mean_even = zscore(mean_even);
    % odd values
    alpha_odd = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'alpha'),'odd'));
    bb_odd  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bb'),'odd'));
    gamma_odd  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'gamma'),'odd'));
    bold_odd  = zscore(ns_mean_by_stimulus(NS, ns_get(NS, 'bold'),'odd'));
    all_mean_data = squeeze(mean(sum(NS.data.ts(:,:,:),2),1));
    mean_odd = zeros(max(NS.trial.condition_num),1);
    for m=1:max(NS.trial.condition_num)+1
        mean_odd(m) = mean(all_mean_data(NS.trial.condition_num==m-1));
    end
    mean_odd = zscore(mean_odd);    
    
    % put in output
    out.bold_vals(k,:) = bold_avg;
    out.bb_vals(k,:) = bb_avg;
    out.gamma_vals(k,:) = gamma_avg;
    out.alpha_vals(k,:) = alpha_avg;
    out.mean_vals(k,:) = mean_avg;
    
    % STATS for different models:
    stats = regstats(bold_even,bb_even); % train 
    bold_pred = stats.beta(1)*ones(size(bold_odd))+stats.beta(2)*bb_odd; % test
    bold_corr(k,1) = corr(bold_pred,bold_odd).^2;
    bold_beta(k,1) = stats.beta(2);

    stats = regstats(bold_even,gamma_even); % train
    bold_pred = stats.beta(1)*ones(size(bold_odd))+stats.beta(2)*gamma_odd; % test
    bold_corr(k,2) = corr(bold_pred,bold_odd).^2;
    bold_beta(k,2) = stats.beta(2);
    
    stats = regstats(bold_even,[bb_even gamma_even]); % train
    bold_pred = stats.beta(1)*ones(size(bold_odd)) + stats.beta(2)*bb_odd  + stats.beta(3)*gamma_odd; % test
    bold_corr(k,3) = corr(bold_pred,bold_odd).^2;
    bold_beta(k,3:4) = stats.beta(2:3);

    stats = regstats(bold_even,alpha_even);% train
    bold_pred = stats.beta(1)*ones(size(bold_odd))+stats.beta(2)*alpha_odd; % test
    bold_corr(k,4) = corr(bold_pred,bold_odd).^2;
    bold_beta(k,5) = stats.beta(2);
    
    stats = regstats(bold_even,[bb_even alpha_even]); % train
    bold_pred = stats.beta(1)*ones(size(bold_odd)) + stats.beta(2)*bb_odd  + stats.beta(3)*alpha_odd; % test
    bold_corr(k,5) = corr(bold_pred,bold_odd).^2;
    bold_beta(k,[6:7]) = stats.beta(2:3);
    
    stats = regstats(bold_even,[gamma_even alpha_even]); %train
    bold_pred = stats.beta(1)*ones(size(bold_odd)) + stats.beta(2)*gamma_odd  + stats.beta(3)*alpha_odd; % test
    bold_corr(k,6) = corr(bold_pred,bold_odd).^2;
    bold_beta(k,[8:9]) = stats.beta(2:3);

    stats = regstats(bold_even,[bb_even gamma_even alpha_even]);% train
    bold_pred = stats.beta(1)*ones(size(bold_odd)) + stats.beta(2)*bb_odd + stats.beta(3)*gamma_odd  + stats.beta(4)*alpha_odd; % test
    bold_corr(k,7) = corr(bold_pred,bold_odd).^2;
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

%% Plot beta and r2 for each model separately
bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};
figure('Position',[0 0 580 200])

for s = 1:length(ns_files)
    subplot(1,2,s),hold on
    for k=1:7
        bar(k,bold_corr(s,k),.9,'FaceColor',bar_colors{k})
    end
    xlim([0 8]),ylim([0 1])
    set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb g','a','bb a','g a','bb g a'})
    title(['Model ' int2str(s) ' R^2'])
end
set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/r_varyBB_sim' int2str(sim_nr)])
% print('-depsc','-r300',['./figures/r_varyBB_sim' int2str(sim_nr)])

figure('Position',[0 0 220 200])
labels_beta={{'bb','',''},{'','g',''},{'bb','g',''},{'','','a'},...
    {'bb','','a'},{'','g','a'},{'bb','g','a'}};
labels_index={[1],[2],[1 2],[3],[1 3],[2 3],[1 2 3]};
bb_g_a_color={[.9 .9 .9],[.6 .6 .6],[.3 .3 .3]};

for s = 1:length(ns_files)
    for k=1:7
        xl_ind=labels_index{k};
        subplot(length(ns_files),7,(s-1)*7+k),hold on 

        temp_beta=bold_beta(s,beta_ind==k);      
        for m=1:length(xl_ind)
            bar(xl_ind(m),temp_beta(m),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
        end

        xlim([.5 3.5]),ylim([-1 1.5])
        set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-1:.2:1.5],'YTickLabel',[])
    end
end
set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/beta_varyBB_sim' int2str(sim_nr)])
% print('-depsc','-r300',['./figures/beta_varyBB_sim' int2str(sim_nr)])

%% correlation across all frequencies

figure('Position',[0 0 400 120])
for s = 1:length(ns_files)
    subplot(1,2,s),hold on
    plot(f,zeros(size(f)),'Color',[.5 .5 .5])
    plot(f,r_all(s,:),'k','LineWidth',2)
    xlabel('Frequency (Hz)'), ylabel('correlation (r)')
    xlim([0 200])
    ylim([-1 1])
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/allfreq_corr_nr' int2str(sim_nr)])
print('-depsc','-r300',['./figures/allfreq_corr_nr' int2str(sim_nr)])


%% plot inputs table
for sim_nr=1:length(NSall)
    NS = NSall{sim_nr};
    plot_colors = [0 0 0; jet(NSall{sim_nr}.params.num_conditions)];

    figure('Position',[0 0 80 280])
    subplot(6,1,1),hold on
    for k=1:length(NS.params.poisson_bb)
        bar(k,NS.params.poisson_baseline+NS.params.poisson_bb(k),'FaceColor',plot_colors(k,:))
    end
    ylim([.4 1])
    set(gca,'YTick',[0:.2:2])
    ylabel('bb level')
    subplot(6,1,2),hold on
    for k=1:length(NS.params.poisson_g)
        bar(k,0+NS.params.poisson_g(k),'FaceColor',plot_colors(k,:))
    end
    ylabel('gamma level')
    ylim([0 1])
    subplot(6,1,3),hold on
    for k=1:length(NS.params.poisson_a)
        bar(k,0+NS.params.poisson_a(k),'FaceColor',plot_colors(k,:))
    end
    ylabel('alpha level')
    ylim([0 .11])

    subplot(6,1,4),hold on
    for k=1:length(NS.params.poisson_bb)
        bar(k,0,'FaceColor',plot_colors(k,:))
    end
    ylabel('bb coh')
    ylim([-0.1 1])
    subplot(6,1,5),hold on
    for k=1:length(NS.params.coherence_g)
        bar(k,NS.params.coherence_g(k),'FaceColor',plot_colors(k,:))
    end
    ylabel('gamma coh')
    ylim([-0.1 1.1])
    subplot(6,1,6),hold on
    for k=1:length(NS.params.poisson_a)
        bar(k,.5,'FaceColor',plot_colors(k,:))
    end
    ylim([-0.1 1])
    ylabel('alpha amp')


    for k=1:6
        subplot(6,1,k)
        box off
        xlim([0 9])
    %     set(gca,'XTick',[1:8])
        set(gca,'XTick',[])
    end
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',['./figures/inputsSim' int2str(sim_nr)])
    print('-depsc','-r300',['./figures/inputsSim' int2str(sim_nr)])

end

%% figure power spectra and BOLD-LFP
%%
for sim_nr=1:length(NSall)
    NS = NSall{sim_nr};

    [bb_avg,bb_std]         = ns_mean_by_stimulus(NS, ns_get(NS, 'bb'));
    [bold_avg,bold_std]     = ns_mean_by_stimulus(NS, ns_get(NS, 'bold'));
    [lfp_avg,lfp_std]       = ns_mean_by_stimulus(NS, ns_get(NS, 'lfp'));
    [gamma_avg,gamma_std]   = ns_mean_by_stimulus(NS, ns_get(NS, 'gamma'));
    [alpha_avg,alpha_std]   = ns_mean_by_stimulus(NS, ns_get(NS, 'alpha'));
    num_conditions          = ns_get(NS, 'num_conditions');
    freq_bb                 = ns_get(NS, 'freq_bb');
    f                       = ns_get(NS, 'f');
    pxx                     = ns_mean_by_stimulus(NS, ns_get(NS, 'power'));
    
    fH=figure('Position',[0 0 900 700]); clf; set(fH, 'Color', 'w')
    fs = [18 12]; % fontsize

    % ---- Plot Spectra for different stimuli -----
    subplot(5,7,[1]), set(gca, 'FontSize', 10);
    plot_colors = [0 0 0; jet(num_conditions)];
    set(gca, 'ColorOrder', plot_colors); hold all
    plot(f(f<max(freq_bb)),pxx(f<max(freq_bb),:), '-', ...
        freq_bb, exp(ns_get(NS, 'power_law')), 'k-', 'LineWidth', 2);
    set(gca, 'XScale', 'log', 'YScale', 'log','XTick',[10 100],'YTick',[10.^-2 10.^0 10.^2])
    xlabel ('Frequency')
    ylabel('Power')
    % xlim([min(freq_bb) max(freq_bb)]);
    xlim([5 max(freq_bb)]);
    ylim([10.^-3 10.^2]);

    % ---- Plot BOLD v ECoG measures ----------------
    num_subplots = 4; % broadband; total LFP; gamma; alpha
    x_data = {bb_avg, lfp_avg, gamma_avg, alpha_avg};
    x_err = {bb_std, lfp_std, gamma_std, alpha_std};
    for ii=1:length(x_err)
        x_err{ii} = x_err{ii}./sqrt(NS.params.num_averages);
    end
    bold_err = bold_std./sqrt(NS.params.num_averages);
    xl     = {'broadband', 'Total LFP power', 'Gamma', 'Alpha'};
    for ii = 1:num_subplots
        this_subplot = 3*(ii-1)+2;
        subplot(num_subplots,3,this_subplot), set(gca, 'FontSize', fs(2)); hold on
        p = polyfit(x_data{ii}, bold_avg,1);
        error_x = [x_data{ii}-x_err{ii} x_data{ii}+x_err{ii}];
        error_y = [bold_avg-bold_err bold_avg+bold_err];
        plot([x_data{ii} x_data{ii}]',error_y','Color',[.5 .5 .5]);
        plot(error_x',[bold_avg bold_avg]','Color',[.5 .5 .5]);
        scatter(x_data{ii}, bold_avg,40,plot_colors(2:end,:)), axis tight square

        hold on; plot(x_data{ii}, polyval(p, x_data{ii}), 'k-', 'LineWidth', 1)
        xlabel(xl{ii},'FontSize',10), ylabel('BOLD','FontSize',10)
        title(sprintf('r = %4.2f', corr(x_data{ii}, bold_avg)));
    end

    % ---- Plot BOLD and ECoG measures as function of simulation inputs -----
    num_subplots = 3; % broadband; total LFP; gamma; alpha
    x_data_name = {'poisson_bb', 'coherence_g', 'poisson_a'};
    xl     = {'Broadband', 'Gamma', 'Alpha'};

    for ii = 1:num_subplots
        this_subplot = 3 * (ii-1)+3;
        subplot(num_subplots,3,this_subplot), set(gca, 'FontSize', fs(2));
        x_data = ns_get(NS, x_data_name{ii});

        [~, inds] = sort(x_data);
        plot(...
            x_data(inds), zscore(bold_avg(inds)), 'go-',...
            x_data(inds), zscore(bb_avg(inds)),  'kd-',...
            x_data(inds), zscore(lfp_avg(inds)),  'bs-',...
            x_data(inds), zscore(gamma_avg(inds)), 'rx-',...
            x_data(inds), zscore(alpha_avg(inds)), 'c*-',...
            'LineWidth', 3)
        xlabel(sprintf('%s (inputs)', xl{ii})), ylabel('response (z-scores)')
        if ii == 1
            legend({'BOLD', 'Broadband',  'LFP power',  'Gamma', 'Alpha'},...
                'Location', 'Best', 'Box', 'off')
        end
    end

    subplot(5,14,3),hold on
    for k=1:8
        bar(k,bold_avg(k),'FaceColor',plot_colors(k,:))
    end
    xlim([0 9]),ylim([0 25])
    set(gca,'XTick',[1:8])
    ylabel('BOLD')

    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',['./figures/powerspectra_Sim' int2str(sim_nr)])
    print('-depsc','-r300',['./figures/powerspectra_Sim' int2str(sim_nr)])

end
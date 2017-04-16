%%
% This script generates pannels for Fig 7 and 8 from Hermes et al:
%
% Purpose: load simulated neural data - time varying membrane potentials -
% that are fit to ECoG data and plot results
%
% DH 2016
%
%% load and plot results for two example electrodes

clear all
sim_nr = 2; % simulation number 2 for paper
elec = 21;%1:22 % examples electrode V1 - 21, V2 - 18

% load the simulation outputs 
%load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_simulation_outputs'],'simulation_outputs')
load(fullfile(BOLD_LFPRootPath, 'data', sprintf('NS_simnr%d_elec%d_simulation_outputs', sim_nr, elec)),'simulation_outputs')

% load the ECoG/fMRI data
load(fullfile(BOLD_LFPRootPath, 'data', 'boldecog_structure_final'))


% load output from the first model (BB - level, G - coh, A - level)
prm_set = 1;
%load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_NS_prmset' int2str(prm_set)],'NS')
load(fullfile(BOLD_LFPRootPath, 'data', sprintf('NS_simnr%d_elec%d_NS_prmset%d', sim_nr, elec,prm_set)),'NS')

data_bb = median(data{elec}.bb_all,2);
data_g = median(data{elec}.gamma_all,2);
data_a = median(data{elec}.alpha_all,2);
data_bold = median(nanmean(data{elec}.allboots,2),3)' * mean(data{elec}.norm);%data{elec}.betas * mean(data{elec}.norm); % to get %signal change
data_bold_ci = data{elec}.se * mean(data{elec}.norm); % ci = 68% confidence interval across bootstraps


%% plot inputs to check
figure('Position',[0 0 200 200])  

plot_colors     = cell2mat(data{elec}.colors(:));
num_conditions  = ns_get(NS, 'num_conditions');

signal_plot = {'bb','g','a'};
for s = 1:length(signal_plot)
    % ---- Plot BB/G/A for different stimuli -----
    subplot(3,2,(2*s-2)+1), set(gca, 'FontSize', 10),hold on
    for ii = 1:num_conditions
        y = getfield(NS.params,['poisson_' signal_plot{s}]);
        if isequal(signal_plot{s},'bb')
            y = y+NS.params.poisson_baseline;
        end
        bar(ii,y(ii),'FaceColor',plot_colors(ii,:))
    end
    ylabel(['poisson ' signal_plot{s}])

    subplot(3,2,(2*s-2)+2), set(gca, 'FontSize', 10),hold on
    for ii = 1:num_conditions
        y = getfield(NS.params,['coherence_' signal_plot{s}]);
        bar(ii,y(ii),'FaceColor',plot_colors(ii,:))
    end
    ylabel(['coherence ' signal_plot{s}])
end
for s = 1:length(signal_plot)
    subplot(3,2,(2*s-2)+1)
    xlim([0 9]),ylim([-.1 1])
    set(gca,'XTick',[1:8])

    subplot(3,2,(2*s-2)+2)
    xlim([0 9]),ylim([-.1 1])
    set(gca,'XTick',[1:8])
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) 'model' int2str(prm_set) '_inputs'])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) 'model' int2str(prm_set) '_inputs'])

%% plot inputs/outputs from a simulation fitted to one ECoG channel

bold_avg        = median(NS.data.bold_bs,2);
bb_avg          = median(NS.data.bb,2);
gamma_avg       = median(NS.data.gamma,2);
alpha_avg       = median(NS.data.alpha,2);
bold_ci         = [quantile(NS.data.bold_bs,.16,2) quantile(NS.data.bold_bs,.84,2)];
bb_ci           = [quantile(NS.data.bb,.16,2) quantile(NS.data.bb,.84,2)];
gamma_ci        = [quantile(NS.data.gamma,.16,2) quantile(NS.data.gamma,.84,2)];
alpha_ci        = [quantile(NS.data.alpha,.16,2) quantile(NS.data.alpha,.84,2)];

num_conditions  = ns_get(NS, 'num_conditions');
f               = ns_get(NS, 'f');
plot_colors     = cell2mat(data{elec}.colors(:));

% INPUTS AND OUTPUTS
figure('Position',[0 0 1100 750])  

% ---- Plot Spectra for different stimuli -----
subplot(4,6,1), set(gca, 'FontSize', 10);
% plot_colors = [0 0 0; jet(num_conditions-1)];
set(gca, 'ColorOrder', plot_colors); hold all

for k = 1:num_conditions
    plot(f,mean(NS.data.lfp_spectra(:,NS.trial.condition_num==k-1),2),...
        'Color',plot_colors(k,:),'LineWidth',2)
end
set(gca, 'XScale', 'log', 'YScale', 'log','XTick',[10 100],'YTick',[10.^-2 10.^0 10.^2])
xlabel ('Frequency'), ylabel('Power')
xlim([5 max(f)]); %ylim([10.^-3 10.^2]);

% ---- Plot BOLD for different stimuli -----
subplot(4,6,2), set(gca, 'FontSize', 10),hold on
for ii = 1:num_conditions
    bar(ii,bold_avg(ii),'FaceColor',plot_colors(ii,:))
end
ylabel('BOLD')

subplot(4,6,3),hold on
for m = 1:length(bold_avg)
    plot(bold_avg(m),data_bold(m),'.','Color',plot_colors(m,:),'MarkerSize',20)
    % plot simulated errorbar (x-axis variance)
    plot(bold_ci(m,:),[data_bold(m) data_bold(m)],'-','Color',[.5 .5 .5])
    % plot data errorbar (y-axis variance)
    plot([bold_avg(m) bold_avg(m)],data_bold_ci(m,:),'-','Color',[.5 .5 .5])
end
r = ns_cod(bold_avg,data_bold', true);
p = polyfit(bold_avg,data_bold',1);
x_line=[min(bold_avg):0.001:max(bold_avg)];
plot(x_line,p(1)*x_line + p(2),'k')
title(['R^2 = ' num2str(r,2)]);
xlim([9.5 15.5]),ylim([-.2 2.2])%ylim([min(data_bold_ci(:))-.2 max(data_bold_ci(:))+.2])
ylabel('measured bold')
xlabel('simulated bold')


% ---- Plot BOLD v ECoG measures ----------------
x_data = {bb_avg, gamma_avg, alpha_avg};
x_err = {bb_ci, gamma_ci, alpha_ci};
xl     = {'broadband', 'gamma', 'alpha'};
for ii = 1:length(x_data)
    subplot(4,6,6+ii), hold on
    p = polyfit(x_data{ii}, bold_avg,1);
    error_x = x_err{ii};
    error_y = bold_ci;
    plot([x_data{ii} x_data{ii}]',error_y','-','Color',[.5 .5 .5]);
    plot(error_x',[bold_avg bold_avg]','-','Color',[.5 .5 .5]);
    scatter(x_data{ii}, bold_avg,40,plot_colors), axis tight square
    hold on; plot(x_data{ii}, polyval(p, x_data{ii}), 'k-', 'LineWidth', 1)
    reg_out = regstats(bold_avg, x_data{ii});
    pred_bold = reg_out.beta(1)+x_data{ii}*reg_out.beta(2:end)';
    xlabel(xl{ii},'FontSize',10), ylabel('BOLD','FontSize',10)
    title(sprintf('R^2 = %4.2f', ns_cod(pred_bold, bold_avg, false)));
end

% ---- Plot BOLD and ECoG measures as function of simulation inputs -----
x_input     = {'poisson_bb', 'poisson_g', 'poisson_a'};
x_output    = {'bb_avg', 'gamma_avg', 'alpha_avg'};
for ii = 1:length(x_data)
    subplot(4,6,12+ii), hold on
    x = getfield(NS.params,x_input{ii});
    y = eval(x_output{ii});
    for k = 1:length(x)
        plot(x(k),y(k),'.','MarkerSize',20,'Color',plot_colors(k,:))
    end
    xlabel(x_input{ii}), ylabel(x_output{ii})
end

x_input     = {'coherence_bb', 'coherence_g', 'coherence_a'};
x_output    = {'bb_avg', 'gamma_avg', 'alpha_avg'};
for ii = 1:length(x_data)
    subplot(4,6,18+ii), hold on
    x = getfield(NS.params,x_input{ii});
    y = eval(x_output{ii});
    for k = 1:length(x)
        plot(x(k),y(k),'.','MarkerSize',20,'Color',plot_colors(k,:))
    end

    xlabel(x_input{ii}), ylabel(x_output{ii})
end

% ---- Plot BOLD correlation with LFP across frequencies -----
subplot(4,3,3), set(gca, 'FontSize', 10),hold on
r_boot=zeros(size(NS.data.lfp_spectra_bs,2),length(f));
for bs=1:size(NS.data.lfp_spectra_bs,2) % number of bootstraps
    fmri_d=NS.data.bold_bs(:,bs);
    ecog_d=squeeze(NS.data.lfp_spectra_bs(:,bs,:));
    r_boot(bs,:)=corr(fmri_d,ecog_d);
end
plot(f,zeros(size(r_boot)),'k','LineWidth',1)
plot(f,median(r_boot,1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])

% REGRESSION MODEL:
bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};

subplot(4,8,14:16),hold on
for k=1:length(NS.stats)
    cross_val_r2 = median(NS.stats(k).stats(:,3),1);
    bar(k,cross_val_r2,.9,'FaceColor',bar_colors{k})
end
xlim([0 8]),ylim([0 1])
set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb g','a','bb a','g a','bb g a'})

labels_beta={{'bb','',''},{'','g',''},{'bb','g',''},{'','','a'},...
    {'bb','','a'},{'','g','a'},{'bb','g','a'}};
labels_index={[1],[2],[1 2],[3],[1 3],[2 3],[1 2 3]};
bb_g_a_color={[.9 .9 .9],[.6 .6 .6],[.3 .3 .3]};

% get minmax beta for plotting max y
min_y = 0;
max_y = 0;
for k=1:length(NS.stats)
    min_y = min([min_y median(NS.stats(k).beta(:,2:end),1)]);
    max_y = max([max_y median(NS.stats(k).beta(:,2:end),1)]);
end

for k=1:length(NS.stats)
    xl_ind=labels_index{k};
    subplot(4,3*length(NS.stats),2*3*length(NS.stats)+2*length(NS.stats)+k),hold on 

    temp_beta=median(NS.stats(k).beta(:,2:end),1);
    for m=1:length(xl_ind)
        bar(xl_ind(m),temp_beta(m),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
    end

%     xlim([.5 3.5]),ylim([min_y-.1 max_y+.1])
    xlim([.5 3.5]),ylim([-.3 .3])
    set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-1:.2:1.5],'YTickLabel',[])
end

% plot LFP simulation versus data
lfp_output    = {'bb', 'gamma', 'alpha'};
data_output    = {'bb_all', 'gamma_all', 'alpha_all'};
for k = 1:3 % bb, g, a
    subplot(4,6,21+k),hold on
    x = mean(getfield(NS.data,lfp_output{k}),2);
    y = mean(getfield(data{elec},data_output{k}),2);
    for m = 1:length(x)
        plot(x(m),y(m),'.','MarkerSize',20,'Color',plot_colors(m,:))
    end
    axis tight
    plot([min(x):.01:max(x)],[min(x):.01:max(x)],'k')
    ylabel(['measured ' lfp_output{k}])
    xlabel(['simulated ' lfp_output{k}])
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) 'model' int2str(prm_set)])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) 'model' int2str(prm_set)])

%
%% plot measured versus predicted LFP and BOLD for all models for one channel

num_conditions  = ns_get(NS, 'num_conditions');
plot_colors = [0 0 0; jet(num_conditions-1)];

figure('Position',[0 0 500 450])
data_bb = median(data{elec}.bb_all,2);
data_g = median(data{elec}.gamma_all,2);
data_a = median(data{elec}.alpha_all,2);
data_bold = data{elec}.betas * mean(data{elec}.norm); % to get %signal change

corr_data_fit = zeros(8,1);
for k=1:8 
    subplot(4,4,k),hold on
    fitted_bold = simulation_outputs(:,k,4);
    for m = 1:length(fitted_bold)
        plot(fitted_bold(m),data_bold(m),'.','Color',data{elec}.colors{m},'MarkerSize',20)
    end
    corr_data_fit(k) = corr(fitted_bold,data_bold');
    p = polyfit(fitted_bold,data_bold',1);
    x_line=[min(fitted_bold):0.001:max(fitted_bold)];
    plot(x_line,p(1)*x_line + p(2),'k')
    title(['R^2 = ' num2str(corr_data_fit(k).^2,2)]);
    xlim([min(min(simulation_outputs(:,:,4)))-.1 max(max(simulation_outputs(:,:,4)))+.1])
    ylim([min(data_bold) max(data_bold)])
    ylabel('measured bold')
    xlabel('simulated bold')
    axis square
end
clear x_line p

for k = 1:8
    subplot(4,4,8+k),hold on
    plot([-1 1],[-1 1],'k')
    plot(simulation_outputs(:,k,1),data_bb,'k.','MarkerSize',10)
    plot(simulation_outputs(:,k,2),data_g,'m.','MarkerSize',10)
    plot(simulation_outputs(:,k,3),data_a,'g.','MarkerSize',10)
%     xlim([-0.1 1]),ylim([-0.1 1])
%     xlim([-1.1 1.1]),ylim([-1.1 1.1])
    axis square
    xlim([-1.5 1.5]),ylim([-1.5 1.5])
    ylabel('measured lfp')
    xlabel('simulated lfp')
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) '_allmodels'])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) '_allmodels'])


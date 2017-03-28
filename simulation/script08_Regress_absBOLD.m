%% now load all fitted electrodes

clear all
sim_nr = 2;
els = 1:1:22;

%%% OUTPUTS:
v_area = NaN(length(els),1); %visual area per electrode
r2_data_fit = NaN(8,length(els)); % R2 between BOLD data and fit for each model:
% DATA: 4:bb-g-a-bold, nr els, up to 10 conditions, 8 simulations
all_data = NaN(4,length(els),10); % ECoG / BOLD data
% SIMULATION: 4:bb-g-a-bold, nr els, up to 10 conditions, 8 simulations
all_simulation = NaN(4,length(els),10,8); % BOLD simulation

% SIMULATION output regression models
all_regressmodels = NaN(length(els),7); % r2 for regression models
all_regressbeta = NaN(length(els),7,4); % betas for regression models

% load the ECoG/fMRI data
load('/Volumes/DoraBigDrive/data/visual/m-files/bold_datalikesimulation/data/boldecog_structure_final.mat');

for l = 1:length(els)
    
    elec = els(l);
       
    % load the simulation outputs 
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_simulation_outputs'],'simulation_outputs')
   
    v_area(l) = data{l}.v_area;

    data_bb = median(data{elec}.bb_all,2);
    data_g = median(data{elec}.gamma_all,2);
    data_a = median(data{elec}.alpha_all,2);
    data_bold = data{elec}.betas * mean(data{elec}.norm);
    
    all_data(1,l,1:length(data_bb))=data_bb;
    all_data(2,l,1:length(data_g))=data_g;
    all_data(3,l,1:length(data_a))=data_a;
    all_data(4,l,1:length(data_bold))=data_bold;
    
    % get simulated ECoG (bb, g, a) and BOLD responses into
    % 'all_simulation' matrix
    for k=1:8 
        % ECoG predictions:
        all_simulation(1,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,1);
        all_simulation(2,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,2);
        all_simulation(3,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,3);
    end

    % load output from the first model (BB - level, G - coh, A - level)
    prm_set = 1;
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_NS_prmset' int2str(prm_set)],'NS')
    % change BOLD output:
    NS = ns_neural2instruments_alternativeBOLD(NS);
    % recalculate bootstrapped BOLD:
    NS = ns_analyse_lfp(NS); %disp(NS.data)
    % recalculate correlations BOLD/LFP:
    NS = ns_summary_statistics(NS);
    
    % put BOLD in for the prem_set 
    all_simulation(4,l,[1:NS.params.num_conditions],prm_set) = median(NS.data.bold_bs,2);
    
    r2_data_fit(prm_set,l) = corr(median(NS.data.bold_bs,2),data_bold').^2;
    
    for k = 1:length(NS.stats)
        % cross validated R3:
        all_regressmodels(l,k) = median(NS.stats(k).stats(:,3));
        % beta values:
        temp_beta = median(NS.stats(k).beta(:,2:end),1);
        all_regressbeta(l,k,1:length(temp_beta)) = temp_beta;
    end
end

%% V1: now plot simulation LFP and BOLD output versus data for all electrodes

cm = lines(length(find(ismember(v_area,1))));

figure('Position',[0 0 700 550])
for k = 1:8
    subplot(4,4,k),hold on
    signal_use = 4; % bold
    x = squeeze(all_simulation(signal_use,ismember(v_area,[1]),:,k));% signal, all electrodes all conditions
    y = squeeze(all_data(signal_use,ismember(v_area,[1]),:)); % signal, all electrodes all conditions
    for l = 1:size(x,1)% regression line for each electrode
        x_el = x(l,:); y_el = y(l,:);
        p = polyfit(x_el(~isnan(x_el)),y_el(~isnan(x_el)),1);
        x_line=[min(x(:)):0.001:max(x(:))];
        plot(x_line,p(1)*x_line + p(2),'Color',cm(l,:))
        plot(x_el,y_el,'.','MarkerSize',10,'Color',cm(l,:))
        clear x_el y_el % housekeeping
    end
    
    % plot each dot with color for condition
%     for m = 1:size(x,2) % plot each value across conditions
%         plot(x(:,m),y(:,m),'.','MarkerSize',10,'Color',data{10}.colors{m})
%     end
    xlim([min(x(:)) max(x(:))]),ylim([min(y(:)) max(y(:))])
    title(['mean R^2 = ' num2str(mean(r2_data_fit(k,ismember(v_area,[1]))),2)])
    xlabel('simulated bold'),ylabel('measured bold')   
end
set(gcf,'PaperPositionMode','auto')

ecog_colors = {'k','r','b'};
for k = 1:8
    subplot(4,4,8+k),hold on
    plot([-1 1],[-1 1],'k')
    for signal_use = 1:3
        x = squeeze(all_simulation(signal_use,ismember(v_area,[1]),:,k));% bb, all electrodes all conditions
        y = squeeze(all_data(signal_use,ismember(v_area,[1]),:)); % bb, all electrodes all conditions
        plot(x,y,'.','MarkerSize',10,'Color',ecog_colors{signal_use})
    end
    axis tight
    xlabel('simulated lfp'),ylabel('measured lfp')   
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/simulatedVSdataLFP_BOLD_allV1electrodes_BOLDabs'])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/simulatedVSdataLFP_BOLD_allV1electrodes_BOLDabs'])

%% V2/V3: now plot simulation LFP and BOLD output versus data for all electrodes
cm = lines(length(find(ismember(v_area,[2 3]))));

figure('Position',[0 0 700 550])
for k = 1:8
    subplot(4,4,k),hold on
    signal_use = 4; % bold
    x = squeeze(all_simulation(signal_use,ismember(v_area,[2 3]),:,k));% signal, all electrodes all conditions
    y = squeeze(all_data(signal_use,ismember(v_area,[2 3]),:)); % signal, all electrodes all conditions
    for l = 1:size(x,1)% regression line for each electrode
        x_el = x(l,:); y_el = y(l,:);
        p = polyfit(x_el(~isnan(x_el)),y_el(~isnan(x_el)),1);
        x_line=[min(x(:)):0.001:max(x(:))];
        plot(x_line,p(1)*x_line + p(2),'Color',cm(l,:))
        plot(x_el,y_el,'.','MarkerSize',10,'Color',cm(l,:))
        clear x_el y_el % housekeeping
    end
    xlim([min(x(:)) max(x(:))]),ylim([min(y(:)) max(y(:))])
    title(['mean R^2 = ' num2str(mean(r2_data_fit(k,ismember(v_area,[2 3]))),2)])
    xlabel('simulated bold'),ylabel('measured bold')   
end
set(gcf,'PaperPositionMode','auto')

ecog_colors = {'k','r','b'};
for k = 1:8
    subplot(4,4,8+k),hold on
    plot([-1 1],[-1 1],'k')
    for signal_use = 1:3
        x = squeeze(all_simulation(signal_use,ismember(v_area,[2 3]),:,k));% bb, all electrodes all conditions
        y = squeeze(all_data(signal_use,ismember(v_area,[2 3]),:)); % bb, all electrodes all conditions
        plot(x,y,'.','MarkerSize',10,'Color',ecog_colors{signal_use})
    end
    axis tight
    xlabel('simulated lfp'),ylabel('measured lfp')   
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/simulatedVSdataLFP_BOLD_allV23electrodes'])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/simulatedVSdataLFP_BOLD_allV23electrodes'])

%% Figure with R2

figure('Position',[0 0 300 400]),hold on
subplot(2,1,1),hold on
for k = [.25 .5 .75]
    plot([0 5],[k k],'Color',[.5 .5 .5])
end
bar(mean(r2_data_fit(1:4,:),2),'FaceColor',[.8 .8 .9])
errorbar([1:4],mean(r2_data_fit(1:4,:),2),std(r2_data_fit(1:4,:),[],2)/sqrt(6),'k.')
ylim([0 1.01])
ylabel('r^2')
xlim([0 5])
set(gca,'YTick',[0:.5:1],'XTick',[1:4],...
    'XTickLabel',{'bLgCaL','bLgCaC','bLgLaL','bLgLaC'})

subplot(2,1,2),hold on
for k =[.25 .5 .75]
    plot([0 5],[k k],'Color',[.5 .5 .5])
end

bar(mean(r2_data_fit(5:8,:),2),'FaceColor',[.8 .8 .9])
errorbar([1:4],mean(r2_data_fit(5:8,:),2),std(r2_data_fit(5:8,:),[],2)/sqrt(6),'k.')
% plot([1:8],corr_data_fit,'k.')
ylim([0 1.01])
ylabel('r^2')
xlim([0 5])
set(gca,'YTick',[0:.5:1],'XTick',[1:4],...
    'XTickLabel',{'bCgCaL','bCgCaC','bCgLaL','bCgLaC'})

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/bestModel'])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/bestModel'])

%%

anova_gr1 = zeros(size(r2_data_fit));
anova_gr1(5:8,:)=1;

anova_gr2 = zeros(size(r2_data_fit));
anova_gr2([3 4 7 8],:)=1;

anova_gr3 = zeros(size(r2_data_fit));
anova_gr3([2 4 6 8],:)=1;

anovan(r2_data_fit(:),{anova_gr1(:) anova_gr2(:) anova_gr3(:)},'model','full')
 


%% R2 plots averaged for V1 and V2 simulations

bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};

figure('Position',[0 0 580 200])
% CROSS-VALIDATED R^2 when taking all boots
subplot(1,2,1),hold on % plot V1
for k=1:size(all_regressmodels,2)
    bar(k,mean(all_regressmodels(v_area==1,k),1),'FaceColor',bar_colors{k})
    % standard error
    mean_resp = mean(all_regressmodels(v_area==1,k),1);
    st_err = std(all_regressmodels(v_area==1,k))./sqrt(sum(ismember(v_area,1)));
    plot([k k],[mean_resp-st_err mean_resp+st_err],'k')
end
clear mean_resp st_err
xlim([0 8]),ylim([0 1])
set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb_g','a','bb_a','g_a','bb_g_a'})
set(gca,'YTick',[0:.2:1])
title('V1 R^2 cross-val')

subplot(1,2,2),hold on % plot V2/V3
for k=1:size(all_regressmodels,2)
    bar(k,mean(all_regressmodels(v_area==2 | v_area==3,k),1),'FaceColor',bar_colors{k})
    % standard error
    mean_resp = mean(all_regressmodels(v_area==2 | v_area==3,k),1);
    st_err = std(all_regressmodels(v_area==2 | v_area==3,k))./sqrt(sum(ismember(v_area,[2 3])));
    plot([k k],[mean_resp-st_err mean_resp+st_err],'k')
end
clear mean_resp st_err
xlim([0 8]),ylim([0 1])
set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb_g','a','bb_a','g_a','bb_g_a'})
set(gca,'YTick',[0:.2:1])
title('V2/V3 R^2 cross-val')

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/r2_plotsV1V23_absBOLD'])
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/r2_plotsV1V23_absBOLD'])


disp(['V1 R^2:' num2str(mean(all_regressmodels(v_area==1,:),1))])
disp(['V2 R^2:' num2str(mean(all_regressmodels(v_area==2 | v_area==3,:),1))])

%% BETA plots averaged for V1 and V2 simulations

labels_beta={{'bb','',''},{'','g',''},{'bb','g',''},{'','','a'},...
    {'bb','','a'},{'','g','a'},{'bb','g','a'}};
labels_index={[1],[2],[1 2],[3],[1 3],[2 3],[1 2 3]};
bb_g_a_color={[.9 .9 .9],[.6 .6 .6],[.3 .3 .3]};

% plot V1
figure('Position',[0 0 450 100])
for k=1:size(all_regressbeta,2)
    xl_ind=labels_index{k};
    subplot(1,size(all_regressbeta,2)*2,k),hold on 
    
    for m=1:3 % nr of betas
        % take the median across the bootstraps for each electrode
        temp_beta=all_regressbeta(v_area==1,k,m);
        if ~isnan(temp_beta(1)) && sum(temp_beta)~=0
            % plot mean across electrodes
            bar(xl_ind(m),mean(temp_beta),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
            % plot 2 x standard error as error bar
            st_err = std(temp_beta)./sqrt(length(temp_beta));
            plot([xl_ind(m) xl_ind(m)],[mean(temp_beta)-st_err mean(temp_beta)+st_err],'k')
            % test for significant difference from zero across electrodes using a t-test
            [~,p]=ttest(temp_beta);
            if p<=0.05
                plot(xl_ind(m),-.2,'r*')
            end
        end
    end
    xlim([.5 3.5]),ylim([-.2 .4])
    set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-0.4:.2:.8],'YTickLabel',[])
end

% plot V2/V3
for k=1:size(all_regressbeta,2)
    xl_ind=labels_index{k};
    subplot(1,size(all_regressbeta,2)*2,size(all_regressbeta,2)+k),hold on 

    for m=1:3 % nr of betas
        % take the median across the bootstraps for each electrode
        temp_beta=all_regressbeta(v_area==2 | v_area==3,k,m);
        if ~isnan(temp_beta(1)) && sum(temp_beta)~=0
            % plot mean across electrodes
            bar(xl_ind(m),mean(temp_beta),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
            % plot 2 x standard error as error bar
            st_err = std(temp_beta)./sqrt(length(temp_beta));
            plot([xl_ind(m) xl_ind(m)],[mean(temp_beta)-st_err mean(temp_beta)+st_err],'k')
            % test for significant difference from zero across electrodes using a t-test
            [~,p]=ttest(temp_beta);
            if p<=0.05
                plot(xl_ind(m),-.2,'r*')
            end
        end
    end   
    xlim([.5 3.5]),ylim([-.2 .4])
    set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-0.4:.2:.8],'YTickLabel',[])
end

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/beta_plotsV1V23_absBOLD'])
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/beta_plotsV1V23_absBOLD'])

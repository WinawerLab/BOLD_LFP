%% now load all fitted electrodes

clear all
sim_nr = 2;
els = 1:1:22;

%%% OUTPUTS:
v_area = NaN(length(els),1); %visual area per electrode
r2_data_fit = NaN(8,length(els)); % R2 between BOLD data and fit for each model:
% DATA: 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations
all_data = NaN(4,length(els),10); % ECoG / BOLD data
% SIMULATION: 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations
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
%         % BOLD baseline subtract and vector length normalize:
%         % subtract baseline:
%         y = simulation_outputs(:,k,4)-simulation_outputs(1,k,4);% subtract baseline
%         y = y/sqrt(sum(y.^2)); % vector length normalize
%         all_simulation(4,l,[1:size(simulation_outputs,1)],k) = y;
        % BOLD, raw from simulation:
        all_simulation(4,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,4);
        
        fitted_bold = simulation_outputs(:,k,4);
        r2_data_fit(k,l) = corr(fitted_bold,data_bold').^2;
    end

    % load output from the first model (BB - level, G - coh, A - level)
    prm_set = 1;
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_NS_prmset' int2str(prm_set)],'NS')
    for k = 1:length(NS.stats)
        % cross validated R3:
        all_regressmodels(l,k) = median(NS.stats(k).stats(:,3));
        % beta values:
        temp_beta = median(NS.stats(k).beta(:,2:end),1);
        all_regressbeta(l,k,1:length(temp_beta)) = temp_beta;
    end
end


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
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/r2_plotsV1V23'])
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/r2_plotsV1V23']) 

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
        if ~isnan(temp_beta)
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
        if ~isnan(temp_beta)
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

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/beta_plotsV1V23'])
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/beta_plotsV1V23'])



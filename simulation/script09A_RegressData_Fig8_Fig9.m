BOLD_LFPRootPath = '/Volumes/DoraBigDrive/github/neural_sim_output/';

load(fullfile(BOLD_LFPRootPath, 'data', 'boldecog_structure_final'));

nr_elec = length(data);
% bounds for quantile plots
ub = .75;
lb = .25;

%% plot BOLD response for one electrode
% el_nr = 14; % 14 = electrode 92
% 
% figure('Position',[0 0 80 150]),hold on
% for k=1:8
%     bar(k,data{el_nr}.betas(k),'FaceColor',data{el_nr}.colors{k})
% end
% xlim([0 9])
% set(gca,'XTick',[1:8])
% set(gcf,'PaperPositionMode','auto')
% title('VL Norm')
% 
% figure('Position',[0 0 80 150]),hold on
% for k=1:8
%     % multiply vector length normalized betas by norm to get an estimate of
%     % % signal change again
%     y = data{el_nr}.betas(k) * mean(data{el_nr}.norm);
%     bar(k,y,'FaceColor',data{el_nr}.colors{k})
% end
% xlim([0 9])
% set(gca,'XTick',[1:8])
% set(gcf,'PaperPositionMode','auto')
% title('%change')
% % print('-dpng','-r300',['./figures/BOLD_data_el' int2str(data{el_nr}.channel)])
% % print('-depsc','-r300',['./figures/BOLD_data_el' int2str(data{el_nr}.channel)])
% clear y

%% Regression analysis ECoG data
% predict ECoGodd/fMRIsubj34 from ECoGeven/fMRIsubj12 and reverse

clear reg_out

% loop regression over electrodes 
v_area=zeros(length(data),1);
for k = 1:2
    reg_out(1,k).stats=NaN(length(data),4); 
    reg_out(2,k).stats=NaN(length(data),4); 
    reg_out(3,k).stats=NaN(length(data),5); 
    reg_out(4,k).stats=NaN(length(data),4); 
    reg_out(5,k).stats=NaN(length(data),5); 
    reg_out(6,k).stats=NaN(length(data),5); 
    reg_out(7,k).stats=NaN(length(data),6); 
    reg_out(8,k).stats=NaN(length(data),4); 
    reg_out(9,k).stats=NaN(length(data),4); 
end

% cross-validated squared Pearson:

r2_crossval_out=NaN(length(data),length(reg_out)); 
% cross-validated coefficient of determination
cod_crossval_out=NaN(length(data),length(reg_out)); 

% fit regression model
for k = 1:nr_elec
    disp(['el ' int2str(k) ' of ' int2str(nr_elec)])
    v_area(k) = data{k}.v_area;
    
    % FIT THE MODEL (training):
    
    % fit model on fmri12, even repeats ECoG; fmri34, odd repeats of ECoG
    fmri_d12 = median(data{k}.allbootsS12,2);
    ecog_bbE = median(data{k}.bb_even,2);
    ecog_gE  = median(data{k}.gamma_even,2);
    ecog_aE  = median(data{k}.alpha_even,2);
   
    fmri_d34 = median(data{k}.allbootsS34,2);
    ecog_bbO = median(data{k}.bb_odd,2);
    ecog_gO  = median(data{k}.gamma_odd,2);
    ecog_aO  = median(data{k}.alpha_odd,2);
    
    % 7 models regressing on different combinations of ECoG measures
    ecog_in{1,1}.data = [ecog_bbE];
    ecog_in{2,1}.data = [ecog_gE];
    ecog_in{3,1}.data = [ecog_bbE ecog_gE];
    ecog_in{4,1}.data = [ecog_aE];
    ecog_in{5,1}.data = [ecog_bbE ecog_aE];
    ecog_in{6,1}.data = [ecog_gE ecog_aE];
    ecog_in{7,1}.data = [ecog_bbE ecog_gE ecog_aE];

    ecog_in{1,2}.data = [ecog_bbO];
    ecog_in{2,2}.data = [ecog_gO];
    ecog_in{3,2}.data = [ecog_bbO ecog_gO];
    ecog_in{4,2}.data = [ecog_aO];
    ecog_in{5,2}.data = [ecog_bbO ecog_aO];
    ecog_in{6,2}.data = [ecog_gO ecog_aO];
    ecog_in{7,2}.data = [ecog_bbO ecog_gO ecog_aO];   
    
    % 1 uniform model
    ecog_in{8,1}.data = [0; ones(size(data{k}.labels,2)-1,1)];    
    ecog_in{8,2}.data = [0; ones(size(data{k}.labels,2)-1,1)];    
    
    % 1 test-retest model
    ecog_in{9,1}.data = [fmri_d34];
    ecog_in{9,2}.data = [fmri_d12];
    
    
    % loop over all models
    for m=1:size(ecog_in,1)
        stats1 = regstats(fmri_d12,ecog_in{m,1}.data); % stats.beta, first one is intercept
        stats2 = regstats(fmri_d34,ecog_in{m,2}.data); % stats.beta, first one is intercept
        
        if ~isnan(stats1.rsquare) % nans if alpha all zeros
            reg_out(m,1).stats(k,1)=stats1.rsquare;
            reg_out(m,1).stats(k,2)=stats1.adjrsquare;
            reg_out(m,1).stats(k,3:2+length(stats1.beta))=stats1.beta; % 1 is the intercept
            reg_out(m,2).stats(k,1)=stats2.rsquare;
            reg_out(m,2).stats(k,2)=stats2.adjrsquare;
            reg_out(m,2).stats(k,3:2+length(stats2.beta))=stats2.beta; % 1 is the intercept
            
        else
            reg_out(m,1).stats(k,1) = 0;
            reg_out(m,1).stats(k,2) = 0;
            reg_out(m,1).stats(k,3:2+length(stats1.beta)) = 0;
            reg_out(m,2).stats(k,1) = 0;
            reg_out(m,2).stats(k,2) = 0;
            reg_out(m,2).stats(k,3:2+length(stats2.beta)) = 0;
                                    
        end
        clear stats1 stats2
    end


    % TEST THE MODEL:

    % CALCULATE PREDICTIONS HERE
    %   We use the regression parameters from the first data set (odd ecog,
    %   fmri S12) to make predictions for the second data set (even ecog,
    %   fmri34) and the opposite - use the regression parameters from the
    %   second data set (even ecog, fmri34) to make predictions for the
    %   first data set

    for m = 1:size(ecog_in,1)
        % Note that for x-validation, we use reg_out(m,2) with ecog_in{m,1}
        % and vice versa
        reg_parms   = reg_out(m,2).stats(k,3:end);
        pred_fmri12 = reg_parms(1)+ecog_in{m,1}.data*reg_parms(2:end)';
                        
        reg_parms   = reg_out(m,1).stats(k,3:end);
        pred_fmri34 = reg_parms(1)+ecog_in{m,2}.data*reg_parms(2:end)';
               
        pred_fmri = [pred_fmri12; pred_fmri34];
        
        if m == 9, pred_fmri = [ecog_in{m,1}.data; ecog_in{m,2}.data]; end
        
        fmri_d    = [fmri_d12; fmri_d34];
        cod_crossval_out(k,m) = ns_cod(pred_fmri,fmri_d); % rescale not necessary, same units
    end
    clear ecog_in
end


%% Shuffled regression analysis

clear reg_outShuff

nr_boot=100; % number of reshuffles 

% loop regression over electrodes and bootstraps
reg_outShuff(1).stats=NaN(length(data),nr_boot,4); 
reg_outShuff(2).stats=NaN(length(data),nr_boot,4); 
reg_outShuff(3).stats=NaN(length(data),nr_boot,5); 
reg_outShuff(4).stats=NaN(length(data),nr_boot,4); 
reg_outShuff(5).stats=NaN(length(data),nr_boot,5); 
reg_outShuff(6).stats=NaN(length(data),nr_boot,5); 
reg_outShuff(7).stats=NaN(length(data),nr_boot,6); 

% for cross-validated R2 and coefficient of determination (takes mean):
r2_crossval_outShuff = NaN(length(data),length(reg_outShuff),nr_boot); 
cod_crossval_outShuff = NaN(length(data),length(reg_outShuff),nr_boot); 

% fit regression model
for k = 1:length(data)
    disp(['el ' int2str(k) ' of ' int2str(length(data))])
    v_area(k) = data{k}.v_area;
    % FIT THE MODEL:
    % fit model on fmri S12, even repeats ECoG
    % train model on fmri S12, odd repeats ECoG
    
    nr_stim = size(data{k}.bb_even,1);
    
    fmri_d34 = median(data{k}.allbootsS34,2);
    fmri_d12 = median(data{k}.allbootsS12,2);

    
    ecog_bbE = median(data{k}.bb_even,2);
    ecog_gE  = median(data{k}.gamma_even,2);
    ecog_aE  = median(data{k}.alpha_even,2);

    ecog_bbO = median(data{k}.bb_odd,2);
    ecog_gO  = median(data{k}.gamma_odd,2);
    ecog_aO  = median(data{k}.alpha_odd,2);
    
    ecog_in{1,1}.data = [ecog_bbE];
    ecog_in{2,1}.data = [ecog_gE];
    ecog_in{3,1}.data = [ecog_bbE ecog_gE];
    ecog_in{4,1}.data = [ecog_aE];
    ecog_in{5,1}.data = [ecog_bbE ecog_aE];
    ecog_in{6,1}.data = [ecog_gE ecog_aE];
    ecog_in{7,1}.data = [ecog_bbE ecog_gE ecog_aE];
    
    ecog_in{1,2}.data = [ecog_bbO];
    ecog_in{2,2}.data = [ecog_gO];
    ecog_in{3,2}.data = [ecog_bbO ecog_gO];
    ecog_in{4,2}.data = [ecog_aO];
    ecog_in{5,2}.data = [ecog_bbO ecog_aO];
    ecog_in{6,2}.data = [ecog_gO ecog_aO];
    ecog_in{7,2}.data = [ecog_bbO ecog_gO ecog_aO];
    
    for bs = 1:nr_boot
        %%%% reshuffle only non-blank conditions
        % fmri_shuffle = [1 randperm(nr_stim-1)+1];
        %           %%%% reshuffle all conditions
        fmri_shuffle = randperm(nr_stim);
                
        fmri_d12_s = median(data{k}.allbootsS12(fmri_shuffle,:),2);
        fmri_d34_s = median(data{k}.allbootsS34(fmri_shuffle,:),2);
       
        for m=1:size(ecog_in,1)
            stats1 = regstats(fmri_d12_s,ecog_in{m,1}.data); % stats.beta, first one is intercept
            stats2 = regstats(fmri_d34_s,ecog_in{m,2}.data); % stats.beta, first one is intercept
            
            if ~isnan(stats1.rsquare) % nans if alpha all zeros
                reg_outShuff(m,1).stats(k,bs,1)=stats1.rsquare;
                reg_outShuff(m,1).stats(k,bs,2)=stats1.adjrsquare;
                reg_outShuff(m,1).stats(k,bs,3:2+length(stats1.beta))=stats1.beta; % 1 is the intercept
                reg_outShuff(m,2).stats(k,bs,1)=stats2.rsquare;
                reg_outShuff(m,2).stats(k,bs,2)=stats2.adjrsquare;
                reg_outShuff(m,2).stats(k,bs,3:2+length(stats2.beta))=stats2.beta; % 1 is the intercept
            else
                reg_outShuff(m,1).stats(k,bs,1) = NaN;
                reg_outShuff(m,1).stats(k,bs,2) = NaN;
                reg_outShuff(m,1).stats(k,bs,3:2+length(stats1.beta)) = NaN;
                reg_outShuff(m,2).stats(k,bs,1) = NaN;
                reg_outShuff(m,2).stats(k,bs,2) = NaN;
                reg_outShuff(m,2).stats(k,bs,3:2+length(stats2.beta)) = NaN;
            end
            clear stats1 stats2
        end        
    end
    
    %%%%% now test the models from reshuffling
 

    for bs = 1:nr_boot
        for m = 1:size(ecog_in,1)
            reg_parms = squeeze(reg_outShuff(m,2).stats(k,bs,3:end));
            pred_fmri12 = reg_parms(1)+ecog_in{m,1}.data*reg_parms(2:end);
            
            reg_parms = squeeze(reg_outShuff(m,1).stats(k,bs,3:end));
            pred_fmri34 = reg_parms(1)+ecog_in{m,2}.data*reg_parms(2:end);
            
            pred_fmri = [pred_fmri12; pred_fmri34];
            fmri_d    = [fmri_d12; fmri_d34];
            r2_crossval_outShuff(k,m,bs) = sign(corr(pred_fmri,fmri_d)) * corr(pred_fmri,fmri_d).^2;
            cod_crossval_outShuff(k,m,bs) = ns_cod(pred_fmri,fmri_d); % rescale not necessary, same units
        end
    end
    clear ecog_in

end


%% PLOT cross-validated COD-R2 across electrodes
% plot reshuffled R2 as an indication of baseline
%
% for the bootstraps:
% 1) median across 100 bootstraps
% 2) variance across electrodes

bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};
box_colors = zeros(length(bar_colors),3);
for ii = 1:length(bar_colors), box_colors(ii,:) = bar_colors{ii}; end

figure('Position',[0 0 580 200])
%plotted_r2 = NaN(length(reg_out),2);

% CROSS-VALIDATED R^2 when taking all boots
for whichAreas = 1:2
    
    subplot(1,2,whichAreas),hold on % plot V1
    if whichAreas == 1
        whichElectrodes = v_area==1;
    else
        whichElectrodes = v_area==2 | v_area==3;
    end
    
    boxplot(cod_crossval_out(whichElectrodes==1,1:7),'Colors', box_colors);%, 'PlotStyle', 'compact');
    for k=1:length(reg_out)-2
        
        % plot R2 from reshuffeling
        plot([k-.4 k+.4],[1 1] * median(median(cod_crossval_outShuff(whichElectrodes,k,:),3),1), ...
            ':','Color',[.5 .5 .5],'LineWidth',2)
        
        % plot R2 from test-retest
        plot([k-.4 k+.4],[1 1] * median(cod_crossval_out(whichElectrodes,9),1), ...
           '-','Color',[.5 .5 .5],'LineWidth',2)
        
    end
    
    clear mean_resp st_err
    xlim([0 8]),ylim([-1 1])
    set(gca,'XTick',[1:7],'XTickLabel',{'bb','g','bb_g','a','bb_a','g_a','bb_g_a'})
    set(gca,'YTick',[-1:.2:1])
    
    if whichAreas == 1, title('V1 R^2 cross-val')
    else, title('V2/V3 R^2 cross-val'); end
    
end

disp(['R^2: ' num2str(median(cod_crossval_out(v_area==1,:),1))]);
disp(['R^2: ' num2str(median(cod_crossval_out(v_area==2 | v_area==3,:),1))]);

fname = fullfile(BOLD_LFPRootPath, 'figures', 'regress_cod_plots_NoRescale');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)

%% make a figure of the betas per electrode

labels_beta={{'bb','',''},{'','g',''},{'bb','g',''},{'','','a'},...
    {'bb','','a'},{'','g','a'},{'bb','g','a'}};
labels_index={[1],[2],[1 2],[3],[1 3],[2 3],[1 2 3]};
bb_g_a_color={[.9 .9 .9],[.6 .6 .6],[.3 .3 .3]};

% plot V1
figure('Position',[0 0 550 100])
for k=1:7 % models to plot
    xl_ind=labels_index{k};
    subplot(1,length(reg_out)*2,k),hold on 
    
    for m=1:size(reg_out(k,1).stats(v_area==1,4:end),2) % nr of betas
        % take the median across the 2 crossvals for each electrode
%         temp_beta=median([reg_out(k,1).stats(v_area==1,3+m) reg_out(k,2).stats(v_area==1,3+m)],2);
        temp_beta=mean([reg_out(k,1).stats(v_area==1,3+m) reg_out(k,2).stats(v_area==1,3+m)],2);
        % plot mean across electrodes
        bar(xl_ind(m),mean(temp_beta),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
        % plot 2 x standard error as error bar
        st_err = std(temp_beta)./sqrt(length(temp_beta));
        plot([xl_ind(m) xl_ind(m)],[mean(temp_beta)-st_err mean(temp_beta)+st_err],'k')
        % test for significant difference from zero across electrodes using a t-test
%         [~,p]=ttest(temp_beta);
%         [p]=signtest(temp_beta);
        p = 2*(.5-abs(.5-mean(bootstrp(10000, @median, temp_beta)>0)));
        if p<=0.05
            plot(xl_ind(m),-.5,'r*')
        end
    end
    xlim([.5 3.5]),ylim([-.5 1])
    set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-0.4:.2:.8],'YTickLabel',[])
end


% plot V2/V3
for k=1:7 % models to plot
    xl_ind=labels_index{k};
    subplot(1,length(reg_out)*2,length(reg_out)+k),hold on 

    for m=1:size(reg_out(k).stats(v_area==2 | v_area==3,4:end),2) % nr of betas
        % take the median across the 2 crossvals for each electrode
%         temp_beta=median([reg_out(k,1).stats(v_area==2 | v_area==3,3+m) reg_out(k,2).stats(v_area==2 | v_area==3,3+m)],2);
        temp_beta=mean([reg_out(k,1).stats(v_area==2 | v_area==3,3+m) reg_out(k,2).stats(v_area==2 | v_area==3,3+m)],2);
        % plot mean across electrodes
        bar(xl_ind(m),mean(temp_beta),.7,'FaceColor',bb_g_a_color{xl_ind(m)})
        % plot 2 x standard error as error bar
        st_err = std(temp_beta)./sqrt(length(temp_beta));
        plot([xl_ind(m) xl_ind(m)],[mean(temp_beta)-st_err mean(temp_beta)+st_err],'k')
        % test for significant difference from zero across electrodes using a t-test
%         [~,p]=ttest(temp_beta);
%         [p]=signtest(temp_beta);
        p = 2*(.5-abs(.5-mean(bootstrp(10000, @median, temp_beta)>0)));
        if p<=0.05
            plot(xl_ind(m),-.5,'r*')
        end
    end   
    xlim([.5 3.5]),ylim([-.5 1])
    set(gca,'XTick',[1:3],'XTickLabel',labels_beta{k},'YTick',[-0.4:.2:.8],'YTickLabel',[])
end

fname = fullfile(BOLD_LFPRootPath, 'figures', 'betas_plot');
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)


%% plot data versus BOLD for all electrodes:

figure('Position',[0 0 300 550])

% choose visual area
v=[2 3];%[2 3]; % can be a number [1] or more [2 3]

% choose ECoG input {1:7} = {bb, g, [bb g], a, [bb a], [g a], [bb g a]};
% only plots pure data here, so only choose models with only
% bb/gamma/alpha, models [1 2 4]
e_in=4;
ecog_names={'bb','g','bb_g','a','bb_a','g_a','bb_g_a'};

spearmanPvals = NaN(1,length(data));
spearmanRhovals = NaN(1,length(data));

v_count=0;
for k=1:length(data)
    disp(['el ' int2str(k) ' of ' int2str(length(data))])
    v_area(k)=data{k}.v_area;
    if ismember(v_area(k),v)
        v_count=v_count+1;
        
        % multiply by the norm to get estimated percent signal change
        fmri_d=median(nanmean(data{k}.allboots,2),3) * mean(data{k}.norm);
        fmri_ci=data{k}.se * mean(data{k}.norm);
        ecog_bbE=median(data{k}.bb_all,2);
        ecog_gE=median(data{k}.gamma_all,2);
        ecog_a=median(data{k}.alpha_all,2);
        
        ecog_in{1}.data=[ecog_bbE];
        ecog_in{1}.ci=[quantile(data{k}.bb_all,lb,2) quantile(data{k}.bb_all,ub,2)];
        ecog_in{2}.data=[ecog_gE];
        ecog_in{2}.ci=[quantile(data{k}.gamma_all,lb,2) quantile(data{k}.gamma_all,ub,2)];
        ecog_in{3}.data=[ecog_bbE ecog_gE];
        ecog_in{4}.data=[ecog_a];
        ecog_in{4}.ci=[quantile(data{k}.alpha_all,lb,2) quantile(data{k}.alpha_all,ub,2)];
        ecog_in{5}.data=[ecog_bbE ecog_a];
        ecog_in{6}.data=[ecog_gE ecog_a];
        ecog_in{7}.data=[ecog_bbE ecog_gE ecog_a];
    
        stats1 = regstats(fmri_d,ecog_in{e_in}.data); % stats.beta, first one is intercept
        % predicted BOLD
        fmri_pred=stats1.beta(1)+ecog_in{e_in}.data*stats1.beta(2:end);
        
        % regression line, sort values just to get them in order
        [x,x_ind]=sort(ecog_in{e_in}.data);
        y=fmri_pred(x_ind);
        
        subplot(5,3,v_count),hold on
        % regression line
        plot(x,y,'k')
        for s=1:length(fmri_d)
            % plot confidence intervals
            plot(ecog_in{e_in}.ci(s,:),[fmri_d(s) fmri_d(s)],'k')
            plot([ecog_in{e_in}.data(s) ecog_in{e_in}.data(s)],fmri_ci(s,:),'k')
            
            % plot data points
            plot(ecog_in{e_in}.data(s),fmri_d(s),'.','Color',data{k}.colors{s},...
                'MarkerSize',20)
        end
        if max(fmri_d)<1
            set(gca,'XTick',[-1.2:.4:1.6],'YTick',[0:.5:5])
        else
            set(gca,'XTick',[-1.2:.4:1.6],'YTick',[0:1:5])
        end

        xlim([min([0; min(ecog_in{e_in}.ci(:))])-.1 max([0; max(ecog_in{e_in}.ci(:))])+.1])
        ylim([min([0; min(fmri_ci(:))])-.1 max([0; max(fmri_ci(:))])+.1]);
        
        title(['R^2 = ' num2str(median(cod_crossval_out(k,e_in,:),3),2)])
        
        % Spearman's rho
        [spearmanRho,spearmanP] = corr(ecog_in{e_in}.data,fmri_d,'Type','Spearman');
        spearmanRhovals(k) = spearmanRho;
        spearmanPvals(k) = spearmanP;
        
    end
end

p = signtest(spearmanRhovals);
disp(['sign test on Spearman rho = ' num2str(p)])

set(gcf,'PaperPositionMode','auto')
fname = fullfile(BOLD_LFPRootPath, 'figures', sprintf('ecogBold_V_%d_%s_nolabels', v(1), ecog_names{e_in}));
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)



%% plot prediction versus BOLD for all electrodes:


figure('Position',[0 0 350 650])

% choose visual area
v=1;%[2 3]; % can be a number [1] or more [2 3]

% choose ECoG input {1:7} = {bb, g, [bb g], a, [bb a], [g a], [bb g a]};
e_in=5;
ecog_names={'bb','g','bb_g','a','bb_a','g_a','bb_g_a'};

v_count=0;
for k=1:length(data)
    disp(['el ' int2str(k) ' of ' int2str(length(data))])
    v_area(k)=data{k}.v_area;
    if ismember(v_area(k),v)
        v_count=v_count+1;
        
        fmri_d=median(nanmean(data{k}.allboots,2),3) * mean(data{k}.norm);
        fmri_ci=data{k}.se * mean(data{k}.norm);       
        ecog_bb=median(data{k}.bb_all,2);
        ecog_g=median(data{k}.gamma_all,2);
        ecog_a=median(data{k}.alpha_all,2);
        
        ecog_in{1}.data=[ecog_bb];
        ecog_in{2}.data=[ecog_g];
        ecog_in{3}.data=[ecog_bb ecog_g];
        ecog_in{4}.data=[ecog_a];
        ecog_in{5}.data=[ecog_bb ecog_a];
        ecog_in{6}.data=[ecog_g ecog_a];
        ecog_in{7}.data=[ecog_bb ecog_g ecog_a];
    
        stats1 = regstats(fmri_d,ecog_in{e_in}.data); % stats.beta, first one is intercept
        % predicted BOLD
        fmri_pred=stats1.beta(1)+ecog_in{e_in}.data*stats1.beta(2:end);
        
        % regression line, sort values just to get them in order
        [b]=regress(fmri_pred,fmri_d);
        x=[min(fmri_pred) max(fmri_pred)];
        y=x*b;
               
        subplot(5,3,v_count),hold on
        
        % regression line
        plot(x,y,'k')
        for s=1:length(fmri_d)
            % plot data points
            plot(fmri_pred(s),fmri_d(s),'.','Color',data{k}.colors{s},...
                'MarkerSize',20)
        end
        
        %  axis tight % we lose one datapoint at the edge...
        xlim([min([0; fmri_pred])-.2 max([0; fmri_pred])+.2]);
        ylim([min([0; fmri_d])-.2 max([0; fmri_d])+.2]);
        
        if max(fmri_pred)<1.5
            set(gca,'Xtick',[0:.5:5],'YTick',[0:.5:5])
        else
            set(gca,'Xtick',[0:1:5],'YTick',[0:1:5])
        end
        
        % cross-validated R2
        title(['R^2 = ' num2str(median(cod_crossval_out(k,e_in,:),3),2)])
    end
end

set(gcf,'PaperPositionMode','auto')
fname = fullfile(BOLD_LFPRootPath, 'figures', sprintf('predBold_V_%d_%s_nolabels', v(1), ecog_names{e_in}));
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)
    

%% plot ECoG versus BOLD for all electrodes:


figure('Position',[0 0 350 650])

% choose visual area
v=[2 3]; % can be a number [1] or more [2 3]

% choose ECoG inemput {1,2 4} = {bb, g, [bb g], a, [bb a], [g a], [bb g a]};
e_in=4;
ecog_names={'bb','g','bb_g','a','bb_a','g_a','bb_g_a'};

v_count=0;
for k=1:length(data)
    disp(['el ' int2str(k) ' of ' int2str(length(data))])
    v_area(k)=data{k}.v_area;
    if ismember(v_area(k),v)
        v_count=v_count+1;
        
        fmri_d=median(nanmean(data{k}.allboots,2),3) * mean(data{k}.norm);
        fmri_ci=data{k}.se * mean(data{k}.norm);       
        ecog_bb=median(data{k}.bb_all,2);
        ecog_g=median(data{k}.gamma_all,2);
        ecog_a=median(data{k}.alpha_all,2);
        
        ecog_in{1}.data=[ecog_bb];
        ecog_in{2}.data=[ecog_g];
        ecog_in{3}.data=[ecog_bb ecog_g];
        ecog_in{4}.data=[ecog_a];
        ecog_in{5}.data=[ecog_bb ecog_a];
        ecog_in{6}.data=[ecog_g ecog_a];
        ecog_in{7}.data=[ecog_bb ecog_g ecog_a];
    
        stats1 = regstats(fmri_d,ecog_in{e_in}.data); % stats.beta, first one is intercept
        
        b(k,1:2) = stats1.beta;
        % regression line, sort values just to get them in order
        x=[min(ecog_in{e_in}.data) max(ecog_in{e_in}.data)];
        y=x*stats1.beta(2)+stats1.beta(1);
               
        subplot(5,3,v_count),hold on
        
        % regression line
        plot(x,y,'k')
        for s=1:length(fmri_d)
            % plot data points
            plot(ecog_in{e_in}.data(s),fmri_d(s),'.','Color',data{k}.colors{s},...
                'MarkerSize',20)
        end
        
        %  axis tight % we lose one datapoint at the edge...
        xlim([min([ecog_in{e_in}.data])-.2 max([ecog_in{e_in}.data])+.2]);
        ylim([min([0; fmri_d])-.2 max([0; fmri_d])+.2]);
        
        if max(ecog_in{e_in}.data)<1.5
            set(gca,'Xtick',[-5:.5:5]);%,'YTick',[0:.5:5])
        else
            set(gca,'Xtick',[-5:1:5]);%,'YTick',[0:1:5])
        end
        
        % cross-validated R2
        title(['R^2 = ' num2str(median(cod_crossval_out(k,e_in,:),3),2)])
    end
end

set(gcf,'PaperPositionMode','auto')
fname = fullfile(BOLD_LFPRootPath, 'figures', sprintf('ecog_V_%d_%s_nolabels', v(1), ecog_names{e_in}));
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)


%% test whether the explained variance is related to the size of the ECoG response

data_var = zeros(length(data),3); % electrodes X bb/g/a
model_r2 = zeros(length(data),4); % electrodes X bb/a/bb&a mode

v_area= zeros(length(data),1);

for k=1:length(data)
    % v_area
    v_area(k) = data{k}.v_area;
    
    %%%% BROADBAND
    % get absolute responses (in log10power - mean)
    abs_response = squeeze(mean(data{k}.bb_all,2));
    % absolute change from baseline
    data_var(k,1) = mean(abs_response(2:end));

    %%%% GAMMA
    % get absolute responses (in log10power - mean)
    abs_response = squeeze(mean(data{k}.gamma_all,2));
    % absolute change from baseline
    data_var(k,2) = mean(abs_response(2:end));

    %%%% ALPHA
    % get absolute responses (in log10power - mean)
    abs_response = squeeze(mean(data{k}.alpha_all,2));
    % absolute change from baseline
    data_var(k,3) = mean(abs_response(2:end));

    %%%% R2 for model BB / A / BB&A
    model_r2(k,:) = median(cod_crossval_out(k,[1 2 4 5],:),3);

end

figure('Position',[0 0 600 300],'Color',[1 1 1])
subplot(2,3,1),hold on
plot(data_var(:,1),model_r2(:,1),'k.','MarkerSize',10)
[r,p]=corr(data_var(:,1),model_r2(:,1));
B = regress(model_r2(:,1),[data_var(:,1) ones(size(data_var(:,1)))]);
x = min(data_var(:,1)):0.01:max(data_var(:,1));
plot(x,B(1)*x+B(2),'k-')
xlabel('size broadband response')
ylabel('R^2 bold/bb')
title(['r = ' num2str(r) ' p = ' num2str(p)])
xlim([0 .5])
% ylim tight
box off

subplot(2,3,2),hold on
plot(data_var(:,2),model_r2(:,2),'k.','MarkerSize',10)
[r,p]=corr(data_var(:,2),model_r2(:,2));
B = regress(model_r2(:,2),[data_var(:,2) ones(size(data_var(:,2)))]);
x = min(data_var(:,2)):0.01:max(data_var(:,2));
plot(x,B(1)*x+B(2),'k-')
xlabel('size gamma response')
ylabel('R^2 bold/g')
title(['r = ' num2str(r) ' p = ' num2str(p)])
xlim([0 .9])
% ylim([-0.5 1])
box off

subplot(2,3,3),hold on
plot(data_var(:,3),model_r2(:,3),'k.','MarkerSize',10)
[r,p]=corr(data_var(:,3),model_r2(:,3));
B = regress(model_r2(:,3),[data_var(:,3) ones(size(data_var(:,3)))]);
x = min(data_var(:,3)):0.01:max(data_var(:,3));
plot(x,B(1)*x+B(2),'k-')
xlabel('size alpha response')
ylabel('R^2 bold/a')
title(['r = ' num2str(r) ' p = ' num2str(p)])
xlim([-.9 0])
% ylim([-0.5 1])
box off

set(gcf,'PaperPositionMode','auto')
fname = fullfile(BOLD_LFPRootPath, 'figures', 'R2_BOLDECoG_ECoGchange');
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)

%% BOLD, bb, gamma, alpha change V1 and V2V3

data_var = zeros(length(data),3); % electrodes X bb/g/a
bold_var = zeros(length(data),1); % electrodes X bb/g/a
v_area= zeros(length(data),1);

for k=1:length(data)
    % v_area
    v_area(k) = data{k}.v_area;
    
    %%%% BOLD
    fmri_d=median(nanmean(data{k}.allboots,2),3) * mean(data{k}.norm);
    % absolute change from baseline
    bold_var(k,1) = mean(fmri_d(2:end));

    %%%% BROADBAND
    % get absolute responses (in log10power - mean)
    abs_response = squeeze(mean(data{k}.bb_all,2));
    % absolute change from baseline
    data_var(k,1) = mean(abs_response(2:end));

    %%%% GAMMA
    % get absolute responses (in log10power - mean)
    abs_response = squeeze(mean(data{k}.gamma_all,2));
    % absolute change from baseline
    data_var(k,2) = mean(abs_response(2:end));

    %%%% ALPHA
    % get absolute responses (in log10power - mean)
    abs_response = squeeze(mean(data{k}.alpha_all,2));
    % absolute change from baseline
    data_var(k,3) = mean(abs_response(2:end));

end

%% ECoG and BOLD plot
figure('Position',[0 0 600 300],'Color',[1 1 1])
signal_in = {'Broadband','Gamma','Alpha'};

for s = 1:3 % signal 
    subplot(1,3,s),hold on
    for v = 1:3
        bar(v,mean(data_var(v_area==v,s)),'w')
        plot(v-.25+[1:size(data_var(v_area==v,s))]/20,data_var(v_area==v,s),'k*')
%         errorbar(v,mean(data_var(v_area==v,s)),std(data_var(v_area==v,s)),'k')
    end
    title(signal_in{s})
    ylabel('Mean change from baseline')
    set(gca,'XTick',[1 2 3],'XTickLabel',{'V1','V2','V3'})
end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/paper_V03/ECoGchangeV123'])
% print('-depsc','-r300',['./figures/paper_V03/ECoGchangeV123'])


figure('Position',[0 0 100 120],'Color',[1 1 1]),hold on
signal_in = {'BOLD'};

for v = 1:3
    bar(v,mean(bold_var(v_area==v)),'w')
    plot(v-.25+[1:size(bold_var(v_area==v))]/20,bold_var(v_area==v),'k*')
%         errorbar(v,mean(data_var(v_area==v,s)),std(data_var(v_area==v,s)),'k')
end
ylabel('Mean change from baseline')
set(gca,'XTick',[1 2 3],'XTickLabel',{'V1','V2','V3'})
[mean(bold_var(v_area==1)); mean(bold_var(v_area==2)); mean(bold_var(v_area==3))]

set(gcf,'PaperPositionMode','auto')
fname = fullfile(BOLD_LFPRootPath, 'figures', 'BOLDchangeV123');
print('-dpng','-r300',fname)
print('-depsc','-r300',fname)


%%
%%
%%
%%
%%
%%
%% NOW CHECKS, CAN BE REMOVED:
%% plot R2 compared to fMRI test-retest and uniform model [0 1 1 1 1 1 1 1]

% plot reshuffled R2 as an indication of baseline
% plot uniform R2 as an indication of baseline

%% load outputs simulation

sim_nr = 2;
%%% OUTPUTS:
r2_data_fit = NaN(8,length(data)); % R2 between BOLD data and fit for each model:
for elec = 1:length(data)
    % get the data
    data_bold = data{elec}.betas * mean(data{elec}.norm);   
    % load the simulation outputs 
    fname = fullfile(BOLD_LFPRootPath, 'data', sprintf('NS_simnr%d_elec%d_simulation_outputs', sim_nr, elec));
    load(fname,'simulation_outputs')   
    for k=1:8 % run across models
        fitted_bold = simulation_outputs(:,k,4);
        r2_data_fit(k,elec) = corr(fitted_bold,data_bold').^2;
    end
end

%% load mean BOLD change
bold_sc = NaN(1,length(data));
for elec = 1:length(data)
    fmri_d=median(nanmean(data{elec}.allboots,2),3) * mean(data{elec}.norm);
    bold_sc(elec) = mean(fmri_d(2:end));
    bold_max(elec) = max(fmri_d);
end

%%
r2_uniform = squeeze(median(cod_crossval_out(:,8,:),3));
r2_testretest = squeeze(median(cod_crossval_out(:,9,:),3));
r2_model5 = squeeze(median(cod_crossval_out(:,5,:),3));
r2_simulation = r2_data_fit(1,:);
el_nrs = [1:length(data)];

figure
subplot(1,2,1),hold on

plot(r2_testretest(v_area==1),'k','LineWidth',2)
plot(r2_uniform(v_area==1),'Color',[.5 .5 .5])
plot(r2_model5(v_area==1),'co','MarkerSize',10)
plot(r2_simulation(v_area==1),'go','MarkerSize',10)
plot(r2_testretest(v_area==1),'k.','MarkerSize',20)
plot(r2_uniform(v_area==1),'.','MarkerSize',20,'Color',[.5 .5 .5])
xlim([0 14]),ylim([0 1])
set(gca,'XTick',[1:length(find(v_area==1))],'XTickLabel',el_nrs(v_area==1))
ylabel('r^2')
title('V1')
legend({'r^2 testretest','r^2 uniform','r^2 regressdata','r^2 simulation'})

subplot(1,2,2),hold on
plot(r2_uniform(v_area==2 | v_area==3),'.','MarkerSize',20,'Color',[.5 .5 .5])
plot(r2_uniform(v_area==2 | v_area==3),'Color',[.5 .5 .5])
plot(r2_testretest(v_area==2 | v_area==3),'k.','MarkerSize',20)
plot(r2_testretest(v_area==2 | v_area==3),'k','LineWidth',2)
plot(r2_model5(v_area==2 | v_area==3),'co','MarkerSize',10)
plot(r2_simulation(v_area==2 | v_area==3),'go','MarkerSize',10)
xlim([0 14]),ylim([0 1])
ylabel('r^2')
title('V2/V3')

set(gca,'XTick',[1:length(find(v_area==2 | v_area==3))],'XTickLabel',el_nrs(v_area==2 | v_area==3))

set(gcf,'PaperPositionMode','auto')
fname = fullfile(BOLD_LFPRootPath, 'figures', sprintf('test_uniform_r2_V1V23'));
print('-painters','-r300','-dpng',fname);
print('-painters','-r300','-depsc',fname);

%% exclude some electrodes with small mean signal change

figure,hold on

t_h = bold_sc'>.7 & r2_testretest>.7;
% t_h = ~(r2_testretest>.7);

plot(r2_testretest(t_h),'k.','MarkerSize',20)
plot(r2_testretest(t_h),'k')
plot(r2_uniform(t_h),'.','MarkerSize',20,'Color',[.5 .5 .5])
plot(r2_uniform(t_h),'Color',[.5 .5 .5])
plot(r2_model5(t_h),'mo','MarkerSize',10)
% plot(r2_simulation(t_h),'go','MarkerSize',10)
ylim([0 1])
el_nrs = [1:length(data)];
set(gca,'XTick',[1:length(find(t_h>0))],'XTickLabel',el_nrs(t_h))



%% plot sorted 

figure('Position',[0 0 500 600])
subplot(2,1,1),hold on
% sorted by simulation r2
[r2_simulation_sorted,sort_els] = sort(r2_simulation,'descend');
v_area_sorted = v_area(sort_els);
plot(r2_testretest(sort_els),'k')
plot(r2_uniform(sort_els),'Color',[.5 .5 .5])
plot(el_nrs(v_area_sorted==1),r2_simulation_sorted(v_area_sorted==1),'r.','MarkerSize',20)
plot(r2_testretest(sort_els),'k.','MarkerSize',20)
plot(r2_uniform(sort_els),'.','MarkerSize',20,'Color',[.5 .5 .5])
% plot(r2_model5(sort_els),'mo','MarkerSize',10)
plot(r2_simulation(sort_els),'go','MarkerSize',10,'LineWidth',3)
ylim([0 1])
el_nrs = [1:length(data)];
set(gca,'XTick',[1:length(sort_els)],'XTickLabel',el_nrs(sort_els))
title(['sorted by r^2 from simulated bold - measured bold'])
legend({'r^2 testretest','r^2 uniform','V1'})
xlabel('electrode nr'),ylabel('r^2')

subplot(2,1,2),hold on
% sorted by ecog-fmri r2
[r2_model5_sorted,sort_els] = sort(r2_model5,'descend');
v_area_sorted = v_area(sort_els);
plot(r2_testretest(sort_els),'k')
plot(r2_uniform(sort_els),'Color',[.5 .5 .5])
plot(el_nrs(v_area_sorted==1),r2_model5_sorted(v_area_sorted==1),'r.','MarkerSize',20)
plot(r2_testretest(sort_els),'k.','MarkerSize',20)
plot(r2_uniform(sort_els),'.','MarkerSize',20,'Color',[.5 .5 .5])
plot(r2_model5(sort_els),'co','MarkerSize',10,'LineWidth',3)
% plot(r2_simulation(sort_els),'go','MarkerSize',10)
ylim([0 1])
el_nrs = [1:length(data)];
set(gca,'XTick',[1:length(sort_els)],'XTickLabel',el_nrs(sort_els))
title(['sorted by r^2 from ECoG-fMRI data model with bb and a'])
xlabel('electrode nr'),ylabel('r^2')

% set(gcf,'PaperPositionMode','auto')
% print('-painters','-r300','-dpng',strcat(['../figures/test_uniform/sorted_r2']));
% print('-painters','-r300','-depsc',strcat(['../figures/test_uniform/sorted_r2']));

%%
figure
hold on
% sorted by simulation r2
[~,sort_els] = sort(bold_sc,'descend');
v_area_sorted = v_area(sort_els);
plot(r2_testretest(sort_els),'k')
plot(r2_uniform(sort_els),'Color',[.5 .5 .5])
plot(r2_model5(sort_els),'co','MarkerSize',10,'LineWidth',3)
plot(r2_simulation(sort_els),'go','MarkerSize',10,'LineWidth',3)
plot(r2_testretest(sort_els),'k.','MarkerSize',20)
plot(r2_uniform(sort_els),'.','MarkerSize',20,'Color',[.5 .5 .5])
ylim([0 1])
el_nrs = [1:length(data)];
set(gca,'XTick',[1:length(sort_els)],'XTickLabel',el_nrs(sort_els))
title(['sorted by bold signal change'])
legend({'r^2 testretest','r^2 uniform','r^2 regress ECoG bb&a','r^2 simulation'})
xlabel('electrode nr'),ylabel('r^2')

set(gcf,'PaperPositionMode','auto')
print('-painters','-r300','-dpng',strcat(['../figures/test_uniform/sorted_by_BOLDsc']));
print('-painters','-r300','-depsc',strcat(['../figures/test_uniform/sorted_by_BOLDsc']));

%% 
clear all
sim_nr = 2;
% load the lookup table from sim_nr:
load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_lookup_table' ],'lookup');
% load one of the calibrations, to get baseline settings correct:
load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_set1' ],'NS')

% load the data
load('/Volumes/DoraBigDrive/data/visual/m-files/bold_datalikesimulation/data/boldecog_structure_final.mat');

%% settings to run the calibrated simulations:

% create table for 8 simulations (orde of this table is g/a/bb)
lookup_combs=[...
    2 3 5
    2 4 5
    1 3 5
    1 4 5
    2 3 6
    2 4 6
    1 3 6
    1 4 6];

nr_conds = 10;
out = script_make_calibration_conds(nr_conds);

poisson_baseline = NS.params.poisson_baseline;
num_neurons = NS.params.num_neurons;
clear NS
%%

tic

for elec = [1:22]%[20 1:4:22] % [1:4:22 20];

nr_conditions = length(data{elec}.stim);

% inputs to simulation
% 8 combinations of bb/gamma/alpha, 3 signals, 2 level/coherence
simulation_inputs = zeros(nr_conditions,8,3,2);

% 8 output values from the simulation (bb,g,a bold)
simulation_outputs = zeros(nr_conditions,8,4);

% lookup input values for 8 conditions, 8 simulations
for cond_nr = 1:nr_conditions

    y_all = [median(data{elec}.bb_even,2) median(data{elec}.gamma_even,2) median(data{elec}.alpha_even,2)];
    y = y_all(cond_nr,:); % 1 electrode, 1 condition
    bb_amp = y(1);
    if bb_amp<0 % can not bb values < 0
        bb_amp = 0;
    end
    gamma_amp = y(2);
    alpha_amp = y(3)+1.5;

    for k = 1:size(lookup_combs,1) % lookup combination
        % broadband
        b = lookup(lookup_combs(k,3)).b; % lookup_combs(k,3) is 5 or 6: vary level/coherence
        if lookup_combs(k,3)==5 % lookup bb level, coherence fixed
            poisson_bb = nthroot(bb_amp./b(1),b(2));
            coherence_bb = out.coherence_bb(1,5);
        elseif lookup_combs(k,3)==6 % lookup bb coherence, level fixed
            poisson_bb = out.poisson_bb(1,6);
            coherence_bb = nthroot(bb_amp./b(1),b(2));
            coherence_bb(coherence_bb>1)=1;
        end

        % gamma
        b = lookup(lookup_combs(k,1)).b; % lookup_combs(k,1) is 1 or 2: vary level/coherence
        if lookup_combs(k,1)==1 % lookup gamma level, coherence fixed
%             poisson_g = nthroot(gamma_amp*(b(2)+bb_amp)./b(1),b(3));
            param1 = (gamma_amp - bb_amp*b(3))*(1+bb_amp);
            if param1<0, param1=0; end
            poisson_g = nthroot(param1./b(1), b(2));
            coherence_g = out.coherence_g(1,1);
        elseif lookup_combs(k,1)==2 % lookup gamma coherence, level fixed
            poisson_g = out.poisson_g(1,2);
            param1 = (gamma_amp - bb_amp*b(3))*(1+bb_amp);
            if param1<0, param1=0; end
            coherence_g = nthroot(param1./b(1), b(2));
            coherence_g(coherence_g>1)=1;
        end

        % alpha
        b = lookup(lookup_combs(k,2)).b; % lookup_combs(k,2) is 3 or 4: vary level/coherence
        if lookup_combs(k,2)==3 % lookup alpha level, coherence fixed
%             poisson_a = nthroot(alpha_amp*(b(2)+bb_amp)./b(1),b(3));
            param1 = (alpha_amp - bb_amp*b(3))*(1+bb_amp);
            if param1<0, param1=0; end
            poisson_a = nthroot(param1./b(1), b(2));
            coherence_a = out.coherence_a(1,3);
        elseif lookup_combs(k,2)==4 % lookup alpha coherence, level fixed
            poisson_a = out.poisson_a(1,4);
            param1 = (alpha_amp - bb_amp*b(3))*(1+bb_amp);
            if param1<0, param1=0; end
            coherence_a = nthroot(param1./b(1), b(2));
            coherence_a(coherence_a>1)=1;
        end
        simulation_inputs(cond_nr,k,1,1)=poisson_bb;
        simulation_inputs(cond_nr,k,1,2)=coherence_bb;
        simulation_inputs(cond_nr,k,2,1)=poisson_g;
        simulation_inputs(cond_nr,k,2,2)=coherence_g;
        simulation_inputs(cond_nr,k,3,1)=poisson_a;
        simulation_inputs(cond_nr,k,3,2)=coherence_a;
    end
end

% run the simulations
for prm_set = 1:size(lookup_combs,1) % 8 combinations
    % Set default parameters
    NS = neural_sim_defaults; %disp(NS.params)
    NS = ns_set(NS, 'save_inputs', 1);
    NS = ns_set(NS, 'num_conditions', size(simulation_inputs,1));
    NS = ns_set(NS, 'num_neurons', num_neurons); 
    NS = ns_set(NS, 'poisson_baseline', poisson_baseline); 
    NS = ns_set(NS, 'poisson_bb',simulation_inputs(:,prm_set,1,1));
    NS = ns_set(NS, 'poisson_g',simulation_inputs(:,prm_set,2,1));
    NS = ns_set(NS, 'poisson_a',simulation_inputs(:,prm_set,3,1));
    NS = ns_set(NS, 'coherence_bb',simulation_inputs(:,prm_set,1,2));
    NS = ns_set(NS, 'coherence_g',simulation_inputs(:,prm_set,2,2));
    NS = ns_set(NS, 'coherence_a',simulation_inputs(:,prm_set,3,2));

    % set bb, gamma and alpha values for all trials
    NS = ns_make_trial_struct(NS);
    % if save_inputs choose trials to save, data can get big if saving a lot
    NS = ns_set(NS, 'trials_save_inputs',[1 length(NS.trial.poisson_bb)]);

    % Simulate. This will produce a time series for each neuron in each trial
    NS = ns_simulate_data(NS); 

    % Convert the neural time series into instrument (LFP/BOLD) measures
    NS = ns_neural2instruments(NS); %disp(NS.data)

    % Get spectra and broadband, gamma, alpha and mean from the LFP
    % note that this needs a baseline, which should be the condition
    % with table 0 in NS.condition_num
    NS = ns_analyse_lfp(NS); %disp(NS.data)

    NS.data.ts = single(NS.data.ts); % to reduce size

    % Compute the correlations between different instrument measures 
    NS = ns_summary_statistics(NS); %disp(NS.stats)
    
    simulation_outputs(:,prm_set,1) = median(NS.data.bb,2); %bb
    simulation_outputs(:,prm_set,2) = median(NS.data.gamma,2); %gamma
    simulation_outputs(:,prm_set,3) = median(NS.data.alpha,2); %alpha
    simulation_outputs(:,prm_set,4) = median(NS.data.bold_bs,2); %bold
    
    save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_NS_prmset' int2str(prm_set)],'NS')
%     load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_NS_prmset' int2str(prm_set)],'NS')
    disp(['DONE elec ' int2str(elec) ' parameter set ' int2str(prm_set)])
    toc
end
% save(['../data/NS_simnr' int2str(sim_nr) '_set' int2str(k) ],'NS')
save(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_simulation_outputs'],'simulation_outputs')
ns_params{k} = NS.params;
toc

end

%% load and plot results for one electrode:

clear all
sim_nr = 2;
for elec = 1:22

% load the simulation outputs 
load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_simulation_outputs'],'simulation_outputs')

% load the ECoG/fMRI data
load('/Volumes/DoraBigDrive/data/visual/m-files/bold_datalikesimulation/data/boldecog_structure_final.mat');

% load output from the first model (BB - level, G - coh, A - level)
prm_set = 1;
load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_NS_prmset' int2str(prm_set)],'NS')


% plot measured versus predicted LFP and BOLD
num_conditions  = ns_get(NS, 'num_conditions');
plot_colors = [0 0 0; jet(num_conditions-1)];

figure('Position',[0 0 700 450])
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
end
clear x_line p

for k = 1:8
    subplot(4,4,8+k),hold on
    plot([-.5 .5],[-.5 .5],'k')
    plot(simulation_outputs(:,k,1),data_bb,'k.','MarkerSize',10)
    plot(simulation_outputs(:,k,2),data_g,'m.','MarkerSize',10)
    plot(simulation_outputs(:,k,3),data_a,'g.','MarkerSize',10)
%     xlim([-0.1 1]),ylim([-0.1 1])
%     xlim([-1.1 1.1]),ylim([-1.1 1.1])
    axis tight
    ylabel('measured lfp')
    xlabel('simulated lfp')
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) '_allmodels'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) '_allmodels'])

% plot inputs/outputs from one specific simulation

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

% INPUTS AND OUTPUTS
figure('Position',[0 0 1100 750])  

% ---- Plot Spectra for different stimuli -----
subplot(4,6,1), set(gca, 'FontSize', 10);
% plot_colors = [0 0 0; jet(num_conditions-1)];
plot_colors = cell2mat(data{elec}.colors(:));
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
fitted_bold = simulation_outputs(:,prm_set,4);
for m = 1:length(fitted_bold)
    plot(fitted_bold(m),data_bold(m),'.','Color',plot_colors(m,:),'MarkerSize',20)
end
r = corr(fitted_bold,data_bold');
p = polyfit(fitted_bold,data_bold',1);
x_line=[min(fitted_bold):0.001:max(fitted_bold)];
plot(x_line,p(1)*x_line + p(2),'k')
title(['R^2 = ' num2str(r.^2,2)]);
axis tight
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
    xlabel(xl{ii},'FontSize',10), ylabel('BOLD','FontSize',10)
    title(sprintf('r = %4.2f', corr(x_data{ii}, bold_avg)));
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

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) 'model' int2str(prm_set)])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) 'model' int2str(prm_set)])

end

%% now load all fitted electrodes

clear all
sim_nr = 2;
els = 1:1:22;
r2_data_fit = NaN(8,length(els));

all_data = NaN(4,length(els),10);
all_simulation = NaN(4,length(els),10,8);
% 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations

for l = 1:length(els)
    
    elec = els(l);

    % load the simulation outputs 
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_simulation_outputs'],'simulation_outputs')

    % load the ECoG/fMRI data
    load('/Volumes/DoraBigDrive/data/visual/m-files/bold_datalikesimulation/data/boldecog_structure_final.mat');

    data_bb = median(data{elec}.bb_all,2);
    data_g = median(data{elec}.gamma_all,2);
    data_a = median(data{elec}.alpha_all,2);
    data_bold = data{elec}.betas * mean(data{elec}.norm);
    
    all_data(1,l,1:length(data_bb))=data_bb;
    all_data(2,l,1:length(data_g))=data_g;
    all_data(3,l,1:length(data_a))=data_a;
    all_data(4,l,1:length(data_bold))=data_bold;
        
    for k=1:8 
        % ECoG:
        all_simulation(1,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,1);
        all_simulation(2,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,2);
        all_simulation(3,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,3);
        % BOLD baseline subtract and vectro length normalize:
        % subtract baseline:
        y = simulation_outputs(:,k,4)-simulation_outputs(1,k,4);% subtract baseline
        y = y/sqrt(sum(y.^2)); % vector length normalize
        all_simulation(4,l,[1:size(simulation_outputs,1)],k) = y;
        % BOLD, raw from simulation:
        all_simulation(4,l,[1:size(simulation_outputs,1)],k) = simulation_outputs(:,k,4);
        
        fitted_bold = simulation_outputs(:,k,4);
        r2_data_fit(k,l) = corr(fitted_bold,data_bold').^2;
    end

end

%% now plot simulation LFP and BOLD output versus data for all electrodes

figure('Position',[0 0 700 500])
for k = 1:8
    subplot(4,4,k),hold on
    signal_use = 4; % bold
    x = squeeze(all_simulation(signal_use,:,:,k));% signal, all electrodes all conditions
    y = squeeze(all_data(signal_use,:,:)); % signal, all electrodes all conditions
    for m = 1:size(x,2)
        plot(x(:,m),y(:,m),'.','MarkerSize',10,'Color',data{10}.colors{m})
    end
    p = polyfit(x(~isnan(x)),y(~isnan(x)),1);
    x_line=[min(x(:)):0.001:max(x(:))];
    plot(x_line,p(1)*x_line + p(2),'k')
    axis tight
    title(['mean R^2 = ' num2str(mean(r2_data_fit(k,:)),2)])
    xlabel('simulated bold'),ylabel('measured bold')   
end
set(gcf,'PaperPositionMode','auto')

ecog_colors = {'k','r','b'};
for k = 1:8
    subplot(4,4,8+k),hold on
    for signal_use = 1:3
        x = squeeze(all_simulation(signal_use,:,:,k));% bb, all electrodes all conditions
        y = squeeze(all_data(signal_use,:,:)); % bb, all electrodes all conditions
        plot(x,y,'.','MarkerSize',10,'Color',ecog_colors{signal_use})
    end
%     xlim([-0.1 1]),ylim([-0.1 1])
%     xlim([-1.1 1.1]),ylim([-1.1 1.1])
    axis tight
    xlabel('simulated lfp'),ylabel('measured lfp')   
end
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/simulatedVSdataLFP_BOLD_allelectrodes'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/simulatedVSdataLFP_BOLD_allelectrodes'])


%%
figure('Position',[0 0 300 400]),hold on
subplot(2,1,1),hold on
for k = [.25 .5 .75]
    plot([0 5],[k k],'Color',[.5 .5 .5])
end
bar(mean(corr_data_fit(1:4,:),2),'FaceColor',[.8 .8 .9])
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

bar(mean(corr_data_fit(5:8,:),2),'FaceColor',[.8 .8 .9])
errorbar([1:4],mean(r2_data_fit(5:8,:),2),std(r2_data_fit(5:8,:),[],2)/sqrt(6),'k.')
% plot([1:8],corr_data_fit,'k.')
ylim([0 1.01])
ylabel('r^2')
xlim([0 5])
set(gca,'YTick',[0:.5:1],'XTick',[1:4],...
    'XTickLabel',{'bCgCaL','bCgCaC','bCgLaL','bCgLaC'})

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/bestModel'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/bestModel'])

%%

anova_gr1 = zeros(size(r2_data_fit));
anova_gr1(5:8,:)=1;

anova_gr2 = zeros(size(r2_data_fit));
anova_gr2([3 4 7 8],:)=1;

anova_gr3 = zeros(size(r2_data_fit));
anova_gr3([2 4 6 8],:)=1;

anovan(corr_data_fit(:),{anova_gr1(:) anova_gr2(:) anova_gr3(:)},'model','full')
 
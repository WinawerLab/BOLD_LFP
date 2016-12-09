%% now load all fitted electrodes

clear all
sim_nr = 2;
els = 1:1:22;

%%% OUTPUTS:
v_area = NaN(length(els),1); %visual area per electrode
% DATA: 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations
all_data = NaN(4,length(els),10); % ECoG / BOLD data
% SIMULATION: 4:bb,g,a,bold, nr els, up to 10 conditions, 8 simulations
all_simulation = NaN(4,length(els),10,8); % BOLD simulation

load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(els(1)) '_NS_prmset' int2str(1)],'NS')
f = ns_get(NS, 'f');
all_r_allfreq = zeros(length(els),length(f));

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
    
    % load output from the first model (BB - level, G - coh, A - level)
    prm_set = 1;
    load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_elec' int2str(elec) '_NS_prmset' int2str(prm_set)],'NS')

    r_boot=zeros(size(NS.data.lfp_spectra_bs,2),length(f));
    for bs=1:size(NS.data.lfp_spectra_bs,2) % number of bootstraps
        fmri_d=NS.data.bold_bs(:,bs);
        ecog_d=squeeze(NS.data.lfp_spectra_bs(:,bs,:));
        r_boot(bs,:)=corr(fmri_d,ecog_d);
    end
    all_r_allfreq(l,:) = median(r_boot,1);   
end

%%

% INPUTS AND OUTPUTS
figure('Position',[0 0 1100 750])  

% ---- Plot BOLD correlation with LFP across frequencies -----
subplot(4,3,3), set(gca, 'FontSize', 10),hold on
plot(f,zeros(size(f)),'k','LineWidth',1)
plot(f,all_r_allfreq(ismember(v_area,1),:),'Color',[.5 .5 .5],'LineWidth',1)
plot(f,all_r_allfreq(21,:),'r','LineWidth',2)
plot(f,mean(all_r_allfreq(ismember(v_area,1),:),1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])

subplot(4,3,6), set(gca, 'FontSize', 10),hold on
plot(f,zeros(size(f)),'k','LineWidth',1)
plot(f,all_r_allfreq(ismember(v_area,[2 3]),:),'Color',[.5 .5 .5],'LineWidth',1)
plot(f,all_r_allfreq(18,:),'r','LineWidth',2)
plot(f,mean(all_r_allfreq(ismember(v_area,[2 3]),:),1),'k','LineWidth',2)
xlabel('Frequency (Hz)'), ylabel('correlation (r)')
xlim([0 200])
ylim([-1 1])

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/Corr_allfreq_model1'])
print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/Corr_allfreq_model1'])


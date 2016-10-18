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
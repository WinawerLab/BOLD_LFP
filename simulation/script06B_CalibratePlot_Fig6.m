%%
% This script generates pannels for Fig 6 from Hermes et al:
%
% Purpose: Load simulated neural data - time varying membrane potentials -
% with a structured set of inputs, such that the relation between the input
% level (amplitude and coherence) of C1, C2 and C3 and the output
% broadband, gamma and alpha power can be described
%
% DH 2016
%
%%
%% Load the lookup table and plot for an example electrode
%%
%% To make figure 6
%%

clear all
sim_nr = 2;
% load the lookup table
load(['/Volumes/DoraBigDrive/github/neural_sim_output/data/NS_simnr' int2str(sim_nr) '_lookup_table' ],'lookup');

% load the data
load('/Volumes/DoraBigDrive/data/visual/m-files/bold_datalikesimulation/data/boldecog_structure_final.mat');
    
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

%% plot lookup table for one electrode

elec = 18; % elec 21 and 18 examples for V1 and V2

nr_conditions = length(data{elec}.stim);
% inputs to simulation
% 8 combinations of bb/gamma/alpha, 3 signals, 2 level/coherence
simulation_inputs = zeros(nr_conditions,8,3,2);

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
        b = lookup(lookup_combs(k,3)).b; % lookup_combs(k,3) is broadband 5 or 6: vary level/coherence
        if lookup_combs(k,3)==5 % lookup bb level, coherence fixed
            poisson_bb = ns_out2in_singleparam(bb_amp,b); % inverse of ns_in2out_singleparam:
            coherence_bb = out.coherence_bb(1,5);
        elseif lookup_combs(k,3)==6 % lookup bb coherence, level fixed
            poisson_bb = out.poisson_bb(1,6);
            coherence_bb = ns_out2in_singleparam(bb_amp,b); % inverse of ns_in2out_singleparam:
            coherence_bb(coherence_bb>1)=1;
            coherence_bb(coherence_bb<0)=0;
        end

        % gamma
        b = lookup(lookup_combs(k,1)).b; % lookup_combs(k,1) is gamma 1 or 2: vary level/coherence
        if lookup_combs(k,1)==1 % lookup gamma level, coherence fixed
            poisson_g = ns_out2in(gamma_amp,b,bb_amp); % inverse of ns_in2out:
            coherence_g = out.coherence_g(1,1);
        elseif lookup_combs(k,1)==2 % lookup gamma coherence, level fixed
            poisson_g = out.poisson_g(1,2);
            coherence_g = ns_out2in(gamma_amp,b,bb_amp);% inverse of ns_in2out:
            coherence_g(coherence_g>1)=1;
            coherence_g(coherence_g<0)=0;
        end

        % alpha
        b = lookup(lookup_combs(k,2)).b; % lookup_combs(k,2) is alpha 3 or 4: vary level/coherence
        if lookup_combs(k,2)==3 % lookup alpha level, coherence fixed
            poisson_a = ns_out2in(alpha_amp,b,bb_amp); % inverse of ns_in2out:
            coherence_a = out.coherence_a(1,3);
        elseif lookup_combs(k,2)==4 % lookup alpha coherence, level fixed
            poisson_a = out.poisson_a(1,4);
            coherence_a = ns_out2in(alpha_amp,b,bb_amp); % inverse of ns_in2out:
            coherence_a(coherence_a>1)=1;
            coherence_a(coherence_a<0)=0;
        end
        simulation_inputs(cond_nr,k,1,1)=poisson_bb;
        simulation_inputs(cond_nr,k,1,2)=coherence_bb;
        simulation_inputs(cond_nr,k,2,1)=poisson_g;
        simulation_inputs(cond_nr,k,2,2)=coherence_g;
        simulation_inputs(cond_nr,k,3,1)=poisson_a;
        simulation_inputs(cond_nr,k,3,2)=coherence_a;
    end
end
clear poisson_bb poisson_g poisson_a coherence_bb coherence_g coherence_a bb_amp gamma_amp alpha_amp

% take the example with bb level, gamma coherence, alpha level (k=1):
%   dimentions: condition, 3 signals, 2 level/coherence
model_nr = 1;
sim_inputs_model1 = squeeze(simulation_inputs(:,model_nr,:,:)); 

% plot broadband poisson in VS LFP broadband out
poisson_bb = sim_inputs_model1(:,1,1);

figure('Position',[0 0 500 150])

subplot(1,5,1:2),hold on

% getted fitted parameters
b = lookup(lookup_combs(model_nr,3)).b; % lookup_combs(1,3) = 5: broadband vary level

% plot fitted function for broadband
y_bbLFP = [-.1:.1:.6];
x_bblevel = ns_out2in_singleparam(y_bbLFP,b); % inverse of ns_in2out_singleparam:
plot(x_bblevel,y_bbLFP,'k')

% plot data points
for cond_nr = 1:length(data{elec}.labels) % stimulus conditions
    % get data broadband value
    bb_avg = median(data{elec}.bb_all(cond_nr,:),2);
    bb_ci = [quantile(data{elec}.bb_all(cond_nr,:),.16,2) quantile(data{elec}.bb_all(cond_nr,:),.84,2)];
    
    % get inputs to create those values
    bb_SimIn = ns_out2in_singleparam(bb_avg,b); % inverse of ns_in2out_singleparam:
    
    % plot lines from x to y
    plot([bb_SimIn max(x_bblevel)],[bb_avg bb_avg],'Color',[.5 1 .5])
    plot([bb_SimIn bb_SimIn],[bb_avg max(y_bbLFP)],'Color',[1 .5 0])

    % plot data points
    plot(bb_SimIn,bb_avg,'.','Color',data{elec}.colors{cond_nr},...
        'MarkerSize',30)
%     plot([bb_SimIn bb_SimIn],bb_ci,'Color',[.5 .5 .5],'LineWidth',2)

end

xlim([min(x_bblevel) max(x_bblevel)]),ylim([min(y_bbLFP) max(y_bbLFP)])
xlabel('C^1 inputs'),ylabel('ECoG broadband')

% C1, C2, C3 bar chart for this electrode
subplot(1,6,4),hold on

%%%% get BB inputs
% measured ECoG bb
bb_avg = median(data{elec}.bb_all,2);
% getted fitted parameters
b = lookup(lookup_combs(model_nr,3)).b; % lookup_combs(1,3) = 5: broadband vary level
x_bblevel = ns_out2in_singleparam(bb_avg,b); % inverse of ns_in2out_singleparam:
% plot data points
for cond_nr = 1:length(data{elec}.labels) % stimulus conditions
    barh(cond_nr,x_bblevel(cond_nr),'FaceColor',data{elec}.colors{cond_nr})
end
xlim([-.1 0.4])

subplot(1,6,5),hold on
%%%% get gamma inputs
% measured ECoG bb and gamma
bb_avg = median(data{elec}.bb_all,2);
gamma_avg = median(data{elec}.gamma_all,2);
% getted fitted parameters
b = lookup(lookup_combs(model_nr,1)).b; % lookup_combs(1,1) = 2: gamma vary coherence
x_gammalevel = ns_out2in(gamma_avg,b,bb_avg); % inverse of ns_in2out_singleparam:
% plot data points
for cond_nr = 1:length(data{elec}.labels) % stimulus conditions
    barh(cond_nr,x_gammalevel(cond_nr),'FaceColor',data{elec}.colors{cond_nr})
end
xlim([-.1 .4])

subplot(1,6,6),hold on
%%%% get alpha inputs
% measured ECoG bb and alpha
bb_avg = median(data{elec}.bb_all,2);
alpha_avg = median(data{elec}.alpha_all,2)+1.5;
% getted fitted parameters
b = lookup(lookup_combs(model_nr,2)).b; % lookup_combs(1,2) = 3: alpha vary level
x_alphalevel = ns_out2in(alpha_avg,b,bb_avg); % inverse of ns_in2out_singleparam:
% plot data points
for cond_nr = 1:length(data{elec}.labels) % stimulus conditions
    barh(cond_nr,x_alphalevel(cond_nr),'FaceColor',data{elec}.colors{cond_nr})
end
xlim([-.1 .6])

for k = 4:6
    subplot(1,6,k),hold on
    set(gca,'YTick',[1:8],'YTickLabel',[])
    ylim([0 9])
    title(['C^' int2str(k-3)])
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) '_bbIn2Out'])
% print('-dpng','-r300',['../figures/sim' int2str(sim_nr) '/Channel' int2str(elec) '_bbIn2Out'])

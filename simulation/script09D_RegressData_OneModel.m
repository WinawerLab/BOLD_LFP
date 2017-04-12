clear all
close all

load(['/Volumes/DoraBigDrive/data/visual/m-files/bold_datalikesimulation/data/boldecog_structure_final.mat'],'data')

%% Put all data in vectors/matrices to analyse in 1 analysis

% calculate matrix size:
nr_rows = 0;
elec_nr = [];
for k = 1:length(data)
    elec_nr(nr_rows+1:nr_rows+length(data{k}.betas)) = k;
    nr_rows = nr_rows + length(data{k}.betas);   
end

% matrices for regression
vals_ecog = NaN(nr_rows,3);
vals_fmri = NaN(nr_rows,1);
v_area = NaN(nr_rows,1);
vect_uniform = NaN(nr_rows,1);

for k = 1:length(data)
    
    v_area(elec_nr==k) = data{k}.v_area;
    
    fmri_d=data{k}.betas;
    ecog_bb=median(data{k}.bb_all,2);
    ecog_g=median(data{k}.gamma_all,2);
    ecog_a=median(data{k}.alpha_all,2);

    % vector length normalize:
    ecog_bb=ecog_bb/sqrt(sum(ecog_bb.^2));
    ecog_g=ecog_g/sqrt(sum(ecog_g.^2));
    ecog_a=ecog_a/sqrt(sum(ecog_a.^2));

    vals_ecog(elec_nr==k,1) = ecog_bb;
    vals_ecog(elec_nr==k,2) = ecog_g;
    vals_ecog(elec_nr==k,3) = ecog_a;
    
    vect_uniform(elec_nr==k) = [0; ones(size(ecog_a,1)-1,1)];
    
    vals_fmri(elec_nr==k) = fmri_d;
end

%% regress with leave-one-out all electrodes together

ecog_models = [1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];

r_vals = zeros(length(data),7);

for k = 1:length(data)
    for m = 1:size(ecog_models,1)
        % training set
        x_matrix = vals_ecog(elec_nr~=k,ecog_models(m,:)>0);
        y_matrix = vals_fmri(elec_nr~=k);

        stats1 = regstats(y_matrix,x_matrix);

        % testing set
        x_matrix = vals_ecog(elec_nr==k,ecog_models(m,:)>0);
        y_matrix = vals_fmri(elec_nr==k,1);

        pred_y = stats1.beta(1) + x_matrix*stats1.beta(2:end);

        r_vals(k,m) = corr(y_matrix,pred_y);
    end
end

%%
figure,hold on
bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};

for k = 1:7
    bar(k,mean(r_vals(:,k).^2),'FaceColor',bar_colors{k})
end
ylim([0 1])
title('all sites')

%% regress with leave-one-out all V1 electrodes

ecog_models = [1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];

vals_ecog_v1 = vals_ecog(v_area==1,:);
vals_fmri_v1 = vals_fmri(v_area==1,:);
elec_nr_v1 = elec_nr(v_area==1);

elec_nrs_incl = unique(elec_nr_v1);

r_vals_v1 = zeros(length(elec_nrs_incl),7);


for k = 1:length(elec_nrs_incl)
    for m = 1:size(ecog_models,1)
        elec = elec_nrs_incl(k);
        % training set
        x_matrix = vals_ecog_v1(elec_nr_v1~=elec,ecog_models(m,:)>0);
        y_matrix = vals_fmri_v1(elec_nr_v1~=elec);

        stats1 = regstats(y_matrix,x_matrix);

        % testing set
        x_matrix = vals_ecog_v1(elec_nr_v1==elec,ecog_models(m,:)>0);
        y_matrix = vals_fmri_v1(elec_nr_v1==elec,1);

        pred_y = stats1.beta(1) + x_matrix*stats1.beta(2:end);

        r_vals_v1(k,m) = corr(y_matrix,pred_y);
    end
end



%% regress with leave-one-out all V2/V3 electrodes

ecog_models = [1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];

vals_ecog_v2 = vals_ecog(v_area==2 | v_area==3,:);
vals_fmri_v2 = vals_fmri(v_area==2 | v_area==3,:);
elec_nr_v2 = elec_nr(v_area==2 | v_area==3);

elec_nrs_incl = unique(elec_nr_v2);

r_vals_v2 = zeros(length(elec_nrs_incl),7);

for k = 1:length(elec_nrs_incl)
    for m = 1:size(ecog_models,1)
        elec = elec_nrs_incl(k);
        % training set
        x_matrix = vals_ecog_v2(elec_nr_v2~=elec,ecog_models(m,:)>0);
        y_matrix = vals_fmri_v2(elec_nr_v2~=elec);

        stats1 = regstats(y_matrix,x_matrix);

        % testing set
        x_matrix = vals_ecog_v2(elec_nr_v2==elec,ecog_models(m,:)>0);
        y_matrix = vals_fmri_v2(elec_nr_v2==elec,1);

        pred_y = stats1.beta(1) + x_matrix*stats1.beta(2:end);

        r_vals_v2(k,m) = corr(y_matrix,pred_y);
    end
end

%% FIGURE
figure
bar_colors={[1 0 0],[1 1 0],[1 .5 0],[0 .2 1],[.5 0 1],[0 .5 0],[.4 .2 .1]};

subplot(1,2,1),hold on
for k = 1:7
    bar(k,mean(r_vals_v1(:,k).^2),'FaceColor',bar_colors{k})
end
ylim([0 1])
title('V1')

subplot(1,2,2),hold on
for k = 1:7
    bar(k,mean(r_vals_v2(:,k).^2),'FaceColor',bar_colors{k})
end
ylim([0 1])
title('V2/V3')

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['../figures/r2_regressAllElectrodesTogether'])
print('-depsc','-r300',['../figures/r2_regressAllElectrodesTogether'])


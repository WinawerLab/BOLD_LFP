function [reg_out,cod_crossval_out,r2_crossval_out] =...
    ns_regress_crossval(fmri_tr,ecog_bb_tr,ecog_g_tr,ecog_a_tr,...
    fmri_te,ecog_bb_te,ecog_g_te,ecog_a_te)

%INPUTS:
% training set
% fmri_tr = median(data{k}.allbootsS12,2);
% ecog_bb_tr = median(data{k}.bb_even,2);
% ecog_g_tr = median(data{k}.gamma_even,2);
% ecog_a_tr = median(data{k}.alpha_even,2);
% 
% % testing set
% fmri_te = median(data{k}.allbootsS34,2);
% ecog_bb_te = median(data{k}.bb_odd,2);
% ecog_g_te = median(data{k}.gamma_odd,2);
% ecog_a_te = median(data{k}.alpha_odd,2);

% outputs
reg_out = [];
r2_crossval_out = NaN(9,1);
cod_crossval_out = NaN(9,1);

ecog_in{1}.data = [ecog_bb_tr];
ecog_in{2}.data = [ecog_g_tr];
ecog_in{3}.data = [ecog_bb_tr ecog_g_tr];
ecog_in{4}.data = [ecog_a_tr];
ecog_in{5}.data = [ecog_bb_tr ecog_a_tr];
ecog_in{6}.data = [ecog_g_tr ecog_a_tr];
ecog_in{7}.data = [ecog_bb_tr ecog_g_tr ecog_a_tr];
% check uniform model
ecog_in{8}.data = [0; ones(size(ecog_bb_tr,1)-1,1)];
% check for fmri_data test-retest
ecog_in{9}.data = [fmri_tr];
        
% estimate m regression models
for m=1:length(ecog_in)
    stats1 = regstats(fmri_tr,ecog_in{m}.data); % stats.beta, first one is intercept
    reg_out(m).stats(1)=stats1.rsquare;
    reg_out(m).stats(2)=stats1.adjrsquare;
    reg_out(m).stats(3:2+length(stats1.beta))=stats1.beta; % 1 is the intercept
    clear stats1
end
clear ecog_in

% TEST THE MODEL:

ecog_in{1}.data = [ecog_bb_te];
ecog_in{2}.data = [ecog_g_te];
ecog_in{3}.data = [ecog_bb_te ecog_g_te];
ecog_in{4}.data = [ecog_a_te];
ecog_in{5}.data = [ecog_bb_te ecog_a_te];
ecog_in{6}.data = [ecog_g_te ecog_a_te];
ecog_in{7}.data = [ecog_bb_te ecog_g_te ecog_a_te];
% check uniform model
ecog_in{8}.data = [0; ones(size(ecog_bb_te,1)-1,1)];
% check fMRI data test-retest
ecog_in{9}.data = fmri_tr;

for m = 1:length(ecog_in)
    reg_parms = reg_out(m).stats(3:end);
    reg_out(m).fmri_pred = reg_parms(1)+ecog_in{m}.data*reg_parms(2:end)';
    x = reg_out(m).fmri_pred;
    y = fmri_te;
    r2_crossval_out(m) = sign(corr(x,y)) * corr(x,y).^2;
    cod_crossval_out(m) = ns_cod(x,y); % rescale not necessary, same units
end


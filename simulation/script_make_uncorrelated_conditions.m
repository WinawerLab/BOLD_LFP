function out = script_make_uncorrelated_conditions(nr_conds)
% Uncorrelate gamma, alpha, from broadband

% to get baselines:
%%%%% vary poisson ranges across sets:
poisson_bb_range =   [.0  .4  ;.0  .15];
poisson_g_range =    [.2  .2  ;.2  .2];
poisson_a_range =    [.2  .5  ;.2  .5];
%%%%% vary coherence ranges across sets:
coherence_bb_range = [.0  .0  ;.0  .0];
coherence_g_range =  [.0  .6  ;.0  .6];
coherence_a_range =  [.75 .75 ;.75 .75];

poisson_bb = NaN(size(poisson_bb_range,1),nr_conds);
poisson_g = NaN(size(poisson_bb_range,1),nr_conds);
poisson_a = NaN(size(poisson_bb_range,1),nr_conds);
coherence_bb = NaN(size(poisson_bb_range,1),nr_conds);
coherence_g = NaN(size(poisson_bb_range,1),nr_conds);
coherence_a = NaN(size(poisson_bb_range,1),nr_conds);

% make the ranges vary:
for k = 1:size(poisson_bb_range,1)
    % Set Poisson rate ranges for all conditions:
    if diff(poisson_bb_range(k,:))==0
        poisson_bb(k,:) = poisson_bb_range(k,1)*ones(1,nr_conds);
    else
        poisson_bb(k,:) = poisson_bb_range(k,1) : diff(poisson_bb_range(k,:))/(nr_conds-1) : poisson_bb_range(k,2);
    end
    
    if diff(poisson_g_range(k,:))==0
        poisson_g(k,:) = poisson_g_range(k,1)*ones(1,nr_conds);
    else
        poisson_g(k,:) = poisson_g_range(k,1) : diff(poisson_g_range(k,:))/(nr_conds-1) : poisson_g_range(k,2);
    end
    if diff(poisson_a_range(k,:))==0
        poisson_a(k,:) = poisson_a_range(k,1)*ones(1,nr_conds);
    else
        poisson_a(k,:) = poisson_a_range(k,1) : diff(poisson_a_range(k,:))/(nr_conds-1) : poisson_a_range(k,2);
    end
    
    % Set coherence ranges for all conditions:
    if diff(coherence_bb_range(k,:))==0
        coherence_bb(k,:) = coherence_bb_range(k,1)*ones(1,nr_conds);
    else
        coherence_bb(k,:) = coherence_bb_range(k,1) : diff(coherence_bb_range(k,:))/(nr_conds-1) : coherence_bb_range(k,2);
    end
    if diff(coherence_g_range(k,:))==0
        coherence_g(k,:) = coherence_g_range(k,1)*ones(1,nr_conds);
    else
        coherence_g(k,:) = coherence_g_range(k,1) : diff(coherence_g_range(k,:))/(nr_conds-1) : coherence_g_range(k,2);
    end
    if diff(coherence_a_range(k,:))==0
        coherence_a(k,:) = coherence_a_range(k,1)*ones(1,nr_conds);
    else
        coherence_a(k,:) = coherence_a_range(k,1) : diff(coherence_a_range(k,:))/(nr_conds-1) : coherence_a_range(k,2);
    end
end

% always have broadband increasing
% and get two additional uncorrelated vectors, here we just assume linearly
% increasing, as set above:
% temporary BB vector: 
temp_bb_vect = 1:1:nr_conds;
% temporary GAMMA vector, to decorrelate from BB 
temp_g_vect = 1:1:nr_conds;
% temporary ALPHA vector, to decorrelate from BB and GAMMA
temp_a_vect = 1:1:nr_conds;

%%%%%%%%%%%%%%%%%%%%%%%
% DECORRELATE temp_gamma_vect from poisson_bb:
b = zeros(100, nr_conds);
c = zeros(nr_conds,1);
for ii = 1:100;
    % define baseline as the minimum gamma
    [~,ind_tmp] = min(temp_g_vect);
    tmp =  [temp_g_vect(ind_tmp) temp_g_vect(setdiff(randperm(nr_conds),ind_tmp(1),'stable'))];
    b(ii,:) = tmp;
    c(ii) = corr(temp_bb_vect', tmp');
end
[~, ind]  = min(abs(c));
ind_g = b(ind,:);
clear tmp ind_tmp ii b c ind

% DECORRELATE temp_gamma_vect from broadband and gamma:
b = zeros(1000, nr_conds);
c = zeros(nr_conds,2); % correlation with bb and gamma
for ii = 1:1000;
    % define baseline as the maximum alpha
    [~,ind_tmp] = max(temp_a_vect);
    tmp =  [temp_a_vect(ind_tmp) temp_a_vect(setdiff(randperm(nr_conds),ind_tmp(1),'stable'))];
    b(ii,:) = tmp;
    c(ii,1) = corr(temp_bb_vect', tmp');
    c(ii,2) = corr(ind_g', tmp');
end
[~, ind]  = min(sum(abs(c),2));
ind_a = b(ind,:);
clear tmp ind_tmp ii b c ind

%%%%%%%%%%%%%%%%%%%%%%%
% now run through all sets and enter decorrelated sequences:
for k = 1:size(poisson_bb_range,1)
    if diff(poisson_g_range(k,:))~=0 % if it should be uncorrelated
        poisson_g(k,:) = poisson_g(k,ind_g);
    end
    if diff(coherence_g_range(k,:))~=0 % if it should be uncorrelated
        coherence_g(k,:) = coherence_g(k,ind_g);
    end
    if diff(poisson_a_range(k,:))~=0 % if it should be uncorrelated
        poisson_a(k,:) = poisson_a(k,ind_a);
    end
    if diff(coherence_a_range(k,:))~=0 % if it should be uncorrelated
        coherence_a(k,:) = coherence_a(k,ind_a);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% create output:
out.poisson_bb = poisson_bb';
out.poisson_g = poisson_g';
out.poisson_a = poisson_a';

out.coherence_bb = coherence_bb';
out.coherence_g = coherence_g';
out.coherence_a = coherence_a';
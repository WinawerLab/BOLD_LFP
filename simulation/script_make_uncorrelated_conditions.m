function out = script_make_uncorrelated_conditions(nr_conds)
% Uncorrelate gamma, alpha, from broadband

%%%%% vary poisson ranges across sets:
poisson_bb_range = [0 .5;0 .2];
poisson_g_range = [.5 .5;.5 .5]; % this should not be 0-0
poisson_a_range = [0 .3;0 .3];
%%%%% vary coherence ranges across sets:
coherence_bb_range = [0 0;0 0];
coherence_g_range = [0 1;0 1];
coherence_a_range = [.5 .5;.5 .5];


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

% uncorrelate poisson_g from poisson_bb in the first set:
k = 1;
b = zeros(100, nr_conds);
c = zeros(nr_conds,1);
for ii = 1:100;
    % define baseline as the minimum poisson_gamma
    [~,ind] = min(poisson_g(k,:));
    tmp =  [poisson_g(k,ind(1)) poisson_g(k, setdiff(randperm(nr_conds),ind(1),'stable'))];
    b(ii,:) = tmp;
    c(ii) = corr(poisson_bb(k,:)', tmp');
end
[~, ind]  = min(abs(c));
poisson_g(k,:) = b(ind,:);

% uncorrelate coherence_g from poisson_bb in the first set:
k = 1;
b = zeros(100, nr_conds);
c = zeros(nr_conds,1);
for ii = 1:100;
    % define baseline as the minimum poisson_gamma
    [~,ind] = min(coherence_g(k,:));
    tmp =  [coherence_g(k,ind(1)) coherence_g(k, setdiff(randperm(nr_conds),ind(1),'stable'))];
    b(ii,:) = tmp;
    c(ii) = corr(poisson_bb(k,:)', tmp');
end
[~, ind]  = min(abs(c));
coherence_g(k,:) = b(ind,:);

% also put in the second set:
k = 2;
coherence_g(k,:) = b(ind,:);

% Decorrelate alpha (in same way as bb/gamma) in the first set:
k = 1;
% Randomize alpha levels across trials
b = zeros(1000, nr_conds);
c = zeros(nr_conds,2); % correlation with bb and gamma
for ii = 1:1000;
    [~,ind] = max(poisson_a(k,:));
    tmp =  [poisson_a(k,ind(1)) poisson_a(k, setdiff(randperm(nr_conds),ind(1),'stable'))];
    b(ii,:) = tmp;
    c(ii,1) = corr(poisson_bb(k,:)', tmp');
    c(ii,2) = corr(coherence_g(k,:)', tmp');
end
[~, ind]  = min(sum(abs(c),2));
poisson_a(k,:) = b(ind,:);

% also put in the second set:
k = 2;
poisson_a(k,:) = b(ind,:);

out.poisson_bb = poisson_bb;
out.poisson_g = poisson_g;
out.poisson_a = poisson_a;

out.coherence_bb = coherence_bb;
out.coherence_g = coherence_g;
out.coherence_a = coherence_a;
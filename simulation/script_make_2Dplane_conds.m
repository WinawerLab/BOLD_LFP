function out = script_make_2Dplane_conds(nr_conds)

%%%%%%%%%%%%%%%%%%%%%%%

% bb coherence:
coh_bb = [0:1/(nr_conds-1):1];

for cal_nr = 1:nr_conds; % loop over bb coherence 
    % vary bb level for each coherence
    inds = ((cal_nr-1)*nr_conds)+[1:nr_conds];
    out.poisson_bb(inds,1) = [0:1/(nr_conds-1):1];
    out.poisson_g(inds,1) = zeros(1,nr_conds);
    out.poisson_a(inds,1) = zeros(1,nr_conds);

    out.coherence_bb(inds,1) = zeros(1,nr_conds)+coh_bb(cal_nr);
    out.coherence_g(inds,1) = zeros(1,nr_conds);
    out.coherence_a(inds,1) = zeros(1,nr_conds);
end


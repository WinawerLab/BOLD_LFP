function out = script_make_calibration_conds(nr_conds)

%%%%%%%%%%%%%%%%%%%%%%%
% create output:
% bb level 0 .5 1, gamma level 0:2
g_level_max = 2;
cal_nr = 1;
out.poisson_bb(:,cal_nr) = [zeros(1,nr_conds) .5+zeros(1,nr_conds) 1+zeros(1,nr_conds)];
out.poisson_g(:,cal_nr) = [0:g_level_max/(nr_conds-1):g_level_max 0:g_level_max/(nr_conds-1):g_level_max 0:g_level_max/(nr_conds-1):g_level_max];
out.poisson_a(:,cal_nr) = zeros(1,30);

out.coherence_bb(:,cal_nr) = zeros(1,30);
out.coherence_g(:,cal_nr) = zeros(1,30);
out.coherence_a(:,cal_nr) = zeros(1,30);

% bb level 0 .5 1, gamma coh 0:1
cal_nr = 2;
out.poisson_bb(:,cal_nr) = [zeros(1,nr_conds) .5+zeros(1,nr_conds) 1+zeros(1,nr_conds)];
out.poisson_g(:,cal_nr) = .2 + zeros(1,30);
out.poisson_a(:,cal_nr) = zeros(1,30);

out.coherence_bb(:,cal_nr) = zeros(1,30);
out.coherence_g(:,cal_nr) = [0:1/(nr_conds-1):1 0:1/(nr_conds-1):1 0:1/(nr_conds-1):1];
out.coherence_a(:,cal_nr) = zeros(1,30);

% bb level 0 .5 1, alpha level 0:5
a_level_max = .5;
cal_nr = 3;
out.poisson_bb(:,cal_nr) = [zeros(1,nr_conds) .5+zeros(1,nr_conds) 1+zeros(1,nr_conds)];
out.poisson_g(:,cal_nr) = zeros(1,30);
out.poisson_a(:,cal_nr) = [a_level_max:-a_level_max/(nr_conds-1):0 a_level_max:-a_level_max/(nr_conds-1):0 a_level_max:-a_level_max/(nr_conds-1):0];

out.coherence_bb(:,cal_nr) = zeros(1,30);
out.coherence_g(:,cal_nr) = zeros(1,30);
out.coherence_a(:,cal_nr) = .75 + zeros(1,30);

% bb level 0 .5 1, alpha coh 0:1
cal_nr = 4;
out.poisson_bb(:,cal_nr) = [zeros(1,nr_conds) .5+zeros(1,nr_conds) 1+zeros(1,nr_conds)];
out.poisson_g(:,cal_nr) = zeros(1,30);
out.poisson_a(:,cal_nr) = .25 + zeros(1,30);

out.coherence_bb(:,cal_nr) = zeros(1,30);
out.coherence_g(:,cal_nr) = zeros(1,30);
out.coherence_a(:,cal_nr) = [1:-1/(nr_conds-1):0 1:-1/(nr_conds-1):0 1:-1/(nr_conds-1):0];

% bb level 0:1
cal_nr = 5;
out.poisson_bb(:,cal_nr) = [0:1/(30-1):1];
out.poisson_g(:,cal_nr) = zeros(1,30);
out.poisson_a(:,cal_nr) = zeros(1,30);

out.coherence_bb(:,cal_nr) = zeros(1,30);
out.coherence_g(:,cal_nr) = zeros(1,30);
out.coherence_a(:,cal_nr) = zeros(1,30);

% bb coherence 0:1
cal_nr = 6;
out.poisson_bb(:,cal_nr) = .2+zeros(1,30);
out.poisson_g(:,cal_nr) = zeros(1,30);
out.poisson_a(:,cal_nr) = zeros(1,30);

out.coherence_bb(:,cal_nr) = [0:1/(30-1):1];
out.coherence_g(:,cal_nr) = zeros(1,30);
out.coherence_a(:,cal_nr) = zeros(1,30);


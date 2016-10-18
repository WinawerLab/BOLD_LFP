function NS = ns_analyse_lfp(NS,varargin)
%
% NS = ns_analyse_lfp(NS)
% 
% the code analyses the simulated data in the same way as the ECoG data

nr_freq         = 200;
num_trials      = ns_get(NS, 'num_trials');
srate           = 1/ns_get(NS,'dt');
alpha_range     = ns_get(NS,'alpha_range');
f_use4fit       = 40:200;
num_conditions  = ns_get(NS, 'num_conditions');

if ~isempty(varargin)
    window_length = varargin{1};
else
    window_length = srate/4;
end

lfp = ns_get(NS, 'lfp');

lfp_spectra = NaN(nr_freq,num_trials);

for k = 1:num_trials
    % get the spectra per trial
    [pxx,f] = ecog_spectra(lfp(:,k),srate,nr_freq,window_length);
    lfp_spectra(:,k) = pxx;    
end

NS = ns_set(NS,'lfp_spectra',lfp_spectra);
NS = ns_set(NS,'f',f);

disp('ns_analyse_lfp has calculated the spectra')

% bootstrap trials
nr_bs = 100;

out.allboots_lfp = NaN(num_conditions,nr_bs,length(f));
out.allbootseven_lfp = NaN(num_conditions,nr_bs,length(f));
out.allbootsodd_lfp = NaN(num_conditions,nr_bs,length(f));
out.allboots_bold = NaN(num_conditions,nr_bs);
out.allbootseven_bold = NaN(num_conditions,nr_bs);
out.allbootsodd_bold = NaN(num_conditions,nr_bs);

for cond = 1:num_conditions
    lfp_spectra_con = lfp_spectra(:,NS.trial.condition_num==cond-1);
    bold_con = NS.data.bold(NS.trial.condition_num==cond-1);
    for bs = 1:nr_bs
        nr_trials_draw = 30;
        % get bootstrap inds for this condition, ALL:
        bs_inds_all = randsample(1:1:size(lfp_spectra_con,2),nr_trials_draw,'true');

        % get bootstrap inds for this condition, EVEN:
        bs_inds_even = randsample(2:2:size(lfp_spectra_con,2),nr_trials_draw,'true');

        % get bootstrap inds for this condition, ODD:
        bs_inds_odd = randsample(1:2:size(lfp_spectra_con,2),nr_trials_draw,'true');

        % bootstrapped spectra
        
        out.allboots_lfp(cond,bs,:) = mean(lfp_spectra_con(:,bs_inds_all),2);
        out.allbootseven_lfp(cond,bs,:) = mean(lfp_spectra_con(:,bs_inds_even),2);
        out.allbootsodd_lfp(cond,bs,:) = mean(lfp_spectra_con(:,bs_inds_odd),2);
        out.allboots_bold(cond,bs) = mean(bold_con(bs_inds_all));
        out.allbootseven_bold(cond,bs) = mean(bold_con(bs_inds_even));
        out.allbootsodd_bold(cond,bs) = mean(bold_con(bs_inds_odd));
    end
end

NS = ns_set(NS, 'lfp_spectra_bs', out.allboots_lfp);
NS = ns_set(NS, 'bold_bs', out.allboots_bold);
NS = ns_set(NS, 'bold_bs_even', out.allbootseven_bold);
NS = ns_set(NS, 'bold_bs_odd', out.allbootsodd_bold);

disp('ns_analyse_lfp has bootstrapped the spectra and BOLD')

%%%%%%% get broadband and gamma per bootstrap

% define the baseline
bs_cond         = 1;
data_base       = squeeze(mean(out.allboots_lfp(bs_cond,:,:),2));

% initialize
bb_weights_all      = zeros(num_conditions,nr_bs);
gamma_weights_all   = zeros(num_conditions,nr_bs);
alpha_weights_all   = zeros(num_conditions,nr_bs);

bb_weights_even     = zeros(num_conditions,nr_bs);
gamma_weights_even  = zeros(num_conditions,nr_bs);
alpha_weights_even  = zeros(num_conditions,nr_bs);

bb_weights_odd      = zeros(num_conditions,nr_bs);
gamma_weights_odd   = zeros(num_conditions,nr_bs);
alpha_weights_odd   = zeros(num_conditions,nr_bs);

for cond = 1:num_conditions
    for bs = 1:nr_bs
        % ALL get bb and gamma
        data_fit                    = squeeze(out.allboots_lfp(cond,bs,:));
        [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
            fit_gammadata(f,f_use4fit,data_base',data_fit');
        bb_weights_all(cond,bs)     = w_pwr;
        gamma_weights_all(cond,bs)  = w_gauss;
        % ALL get alpha
        spectral_change             = log10(data_fit) - log10(data_base);
        alpha_weights_all(cond,bs)  = mean(spectral_change(f>=alpha_range(1) & f<=alpha_range(2)));
        
        % EVEN get bb and gamma
        data_fit                    = squeeze(out.allbootseven_lfp(cond,bs,:));
        [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
            fit_gammadata(f,f_use4fit,data_base',data_fit');
        bb_weights_even(cond,bs)    = w_pwr;
        gamma_weights_even(cond,bs) = w_gauss;
        % EVEN get alpha
        spectral_change             = log10(data_fit) - log10(data_base);
        alpha_weights_even(cond,bs) = mean(spectral_change(f>=alpha_range(1) & f<=alpha_range(2)));

        % ODD get bb and gamma
        data_fit                    = squeeze(out.allbootsodd_lfp(cond,bs,:));
        [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
            fit_gammadata(f,f_use4fit,data_base',data_fit');
        bb_weights_odd(cond,bs)     = w_pwr;
        gamma_weights_odd(cond,bs)  = w_gauss;
        % ODD get alpha
        spectral_change             = log10(data_fit) - log10(data_base);
        alpha_weights_odd(cond,bs)  = mean(spectral_change(f>=alpha_range(1) & f<=alpha_range(2)));
        clear spectral_change data_fit
    end
end

disp('ns_analyse_lfp has extraced bb, gamma, alpha measured from bootstraps')

% adjust broadband with baseline offset:
bb_baseline     = median(bb_weights_all(1,:),2);
bb_weights_all  = bb_weights_all - bb_baseline;
bb_weights_even = bb_weights_even - bb_baseline;
bb_weights_odd  = bb_weights_odd - bb_baseline;
    
% set in structure
NS = ns_set(NS, 'bb', bb_weights_all);
NS = ns_set(NS, 'gamma', gamma_weights_all);
NS = ns_set(NS, 'alpha', alpha_weights_all);

NS = ns_set(NS, 'bb_even', bb_weights_even);
NS = ns_set(NS, 'gamma_even', gamma_weights_even);
NS = ns_set(NS, 'alpha_even', alpha_weights_even);

NS = ns_set(NS, 'bb_odd', bb_weights_odd);
NS = ns_set(NS, 'gamma_odd', gamma_weights_odd);
NS = ns_set(NS, 'alpha_odd', alpha_weights_odd);

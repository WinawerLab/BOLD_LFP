coherence_a     = [0:.01:1];
calculated_cov  = NaN(size(coherence_a));
calculated_corr = NaN(size(coherence_a));

for ii = 1:length(coherence_a)

    % set the synchrony for all neurons, no asynchrony if highly
    % coherent, asynchrony varying across alpha-range if less coherence
    onset_asynchrony = randi([0 round(mean(wait_time)*(1-coherence_a(ii)))],[1 num_neurons]);
    alpha_pulses_short = zeros(length(t), num_neurons);
    for k=1:num_neurons
        alpha_pulses_short(:,k) = alpha_pulses(onset_asynchrony(k)+1:onset_asynchrony(k)+length(t),k);
    end
    h = exp(-(t-.075).^2/(2*.02^2));
    alpha_inputs = -conv2(alpha_pulses_short, h', 'full');
    [~,max_ind] = max(h);
    % get the peak of the response at the time of the pulse: 
    alpha_inputs = alpha_inputs(max_ind:max_ind+length(t)-1,:);

    % save average covariance in the alpha inputs:
    if poisson_a(ii)>0
        alpha_cov = corrcoef(alpha_inputs); % covariance matrix normalized
        alpha_cov(tril(ones(size(cov(alpha_inputs))))>0) = NaN;% set lower triangle to NaN
        calculated_corr(ii) = nanmean(alpha_cov(:));

        alpha_cov = cov(alpha_inputs); % covariance matrix normalized
        alpha_cov(tril(ones(size(cov(alpha_inputs))))>0) = NaN;% set lower triangle to NaN
        calculated_cov(ii) = nanmean(alpha_cov(:));
    else % zero if there is no alpha input
        calculated_corr(ii) = 0;
        calculated_cov(ii) = 0;
    end
end

%%
figure
subplot(2,1,1)
plot(coherence_a,calculated_corr)
xlabel('input coherence')
ylabel('alpha input mean corrcoef')
subplot(2,1,2)
plot(coherence_a,calculated_cov)
xlabel('input coherence')
ylabel('alpha input mean cov')

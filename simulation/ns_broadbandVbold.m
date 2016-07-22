n = 200;   % number of neurons
t = 30;   % number of trials
T = 1000;  % number of time points


v = 4+rand(1,t);         % variance for each trial
x = randn(T * n, t);     % data

x = bsxfun(@times, x,v); 
x = reshape(x, [T, n, t]);
x = cumsum(x);

pow = abs(fft(sum(x,2))).^2;        % sum over neurons, compute power
pow = squeeze(pow(1:T/2,:,:));
bb  = exp(mean(log(pow)))';         % mean across frequencies (log normalized)

BOLD = squeeze(sum(var(x),2));      % sum of variance
%LFP  = squeeze(var(sum(x,2)));
LFP = bb;

tbl = table(LFP, BOLD, 'VariableNames',{'BroadbandLFP','BOLD'});
lm1 = fitlm(tbl);
tbl = table(v'.^2, BOLD, 'VariableNames',{'InputVariance','BOLD'});
lm2 = fitlm(tbl);
tbl = table(v'.^2, LFP, 'VariableNames',{'InputVariance','BroadbandLFP'});
lm3 = fitlm(tbl);
figure(1), clf; 

subplot(3,1,1)
scatter(LFP, BOLD), hold on, plot(lm1);
title(sprintf('%4.2f', 100*lm1.Rsquared.Ordinary))

subplot(3,1,2)
scatter(v.^2, BOLD), hold on, plot(lm2);
title(sprintf('%4.2f', 100*lm2.Rsquared.Ordinary))


subplot(3,1,3)
scatter(v.^2, LFP), hold on, plot(lm3);
title(sprintf('%4.2f', 100*lm3.Rsquared.Ordinary))

return
%%
n = 300;
d = zeros(1,n);
o = zeros(1,n);
for ii = 1:n
    x = randn(n);

    c = cov(x');
    d(ii) = sum(diag(c));
    o(ii) = sum(c(:)) - d(ii);
end
figure(2); clf
scatter(d,o)
axis equal
xlabel('Variance')
ylabel('Covariance')





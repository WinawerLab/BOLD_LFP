x = [0:.1:1];
y = 1./(2+5*x);

figure,
plot(x,y)

%%

mu              = zeros(1,num_broadband); % if you add offset here it would get filtered out
sigma           = eye(num_broadband) + (1-eye(num_broadband))* ns_get(NS, 'alpha_coh');
alpha_inputs    = mvnrnd(mu,sigma,length(t));
alpha_inputs    = poisson_rate_a(ii)*filter(alpha_filter, alpha_inputs);        
alpha_envelope  = abs(hilbert(alpha_inputs));
% I am not sure that this is necessary:
% invert the alpha envelope, such that alpha correlates negatively
% with BOLD:
alpha_envelope  = 1./(2+5*alpha_envelope);
alpha_envelope  = filtfilt(low_bf_b, low_bf_a, alpha_envelope); % lowpass filter the envelope
alpha_inputs    = alpha_inputs + alpha_envelope;


%%


num_broadband   = ns_get(NS,'num_broadband');
srate           = 1/(ns_get(NS,'dt'));
t               = ns_get(NS, 't');

d1 = designfilt('bandpassiir','FilterOrder',10, ...
    'HalfPowerFrequency1',8,'HalfPowerFrequency2',15,...
    'SampleRate',srate,'DesignMethod','butter');

d2 = designfilt('lowpassiir','FilterOrder',10, ...
    'HalfPowerFrequency',5,...
    'SampleRate',srate,'DesignMethod','butter');

mu              = zeros(1,num_broadband); % if you add offset here it would get filtered out
sigma           = eye(num_broadband) + (1-eye(num_broadband))* ns_get(NS, 'alpha_coh');
alpha_inputs    = mvnrnd(mu,sigma,length(t));

figure('Position',[0 0 700 200])

rates_in = [0:.2:1];
for k = 1:length(rates_in)
poisson_rate_a  = rates_in(k);

alpha_inputs1    = poisson_rate_a*filtfilt(d1, alpha_inputs);

alpha_envelope  = abs(hilbert(alpha_inputs1));
alpha_envelope1 = .05./(1+25*alpha_envelope);

alpha_envelope2 = filtfilt(d2, alpha_envelope1);

alpha_inputs2    = alpha_inputs1 + alpha_envelope2;

subplot(2,length(rates_in),k),hold on
plot(alpha_inputs1,'Color',[.5 .5 .5])

plot(alpha_envelope,'r')
plot(alpha_envelope1,'c')

plot(alpha_envelope2,'b')
plot(alpha_inputs2,'k')
ylim([-.5 1])

subplot(2,length(rates_in),length(rates_in)+k),hold on
plot(zeros(size(alpha_inputs2)),':','Color',[.5 .5 .5])
plot(alpha_inputs2,'k')
ylim([-.5 1])

end
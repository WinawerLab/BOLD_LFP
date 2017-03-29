%%
% This script generates pannels for Fig 1 from Hermes et al:
%
% This script generates two example signals with two signusoids that are in
% phase or out of phase and shows the different in the power of the sum and
% the sum of the power of these two signals
%
% DH 2016


%% Simmulate some signals 
srate = 1000;
t = 1/srate:1/srate:.5;
sin_freq = 10;

% signal  = sinus + noise
signal1 = sin(sin_freq*2*pi*t)+1/10*randn(size(t));

% signal 2 is in phase
signal2 = sin(sin_freq*2*pi*t)+1/10*randn(size(t));

% signal 3 is 180 out of phase
signal3 = sin(pi+sin_freq*2*pi*t)+1/10*randn(size(t));


%%
%% Make example for figure 1
%%

figure('Position',[0 0 400 250],'Color',[1 1 1])

% signal 1 and 2
subplot(4,4,1), hold on
plot(t,signal1,'k')
title('x_1')

subplot(4,4,5), hold on
plot(t,signal2,'k')
title('x_2')

subplot(4,4,9), hold on
plot(t,(signal1 + signal2).^2,'k')
title('(x_1 + x_2).^2')
ylim([0 5])

subplot(4,4,10), hold on
plot(t,signal1.^2 + signal2.^2,'k')
title('x_1.^2 + x_2.^2')
ylim([0 5])


% signal 1 and 2
subplot(4,4,3), hold on
plot(t,signal1,'k')
title('x_1')

subplot(4,4,7), hold on
plot(t,signal3,'k')
title('x_2')

subplot(4,4,11), hold on
plot(t,(signal1 + signal3).^2,'k')
title('(x_1 + x_2)^2')
ylim([0 5])

subplot(4,4,12), hold on
plot(t,signal1.^2 + signal3.^2,'k')
title('x_1^2 + x_2^2')
ylim([0 5])

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['../figures/ns_ConceptMethodsFig01'])
% print('-depsc','-r300',['../figures/ns_ConceptMethodsFig01'])

disp(['power(X1)+power(X2))= ' num2str(mean(signal1.^2) + mean(signal2.^2))])
disp(['power(X1)+power(X3))= ' num2str(mean(signal1.^2) + mean(signal3.^2))])
disp(['power(X1+X2)= ' num2str(mean((signal1 + signal2).^2))])
disp(['power(X1+X3)= ' num2str(mean((signal1 + signal3).^2))])

 

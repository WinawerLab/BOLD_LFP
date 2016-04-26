

%% 
srate = 1000;
t = 1/srate:1/srate:1;
sin_freq = 10;

% signal  = sinus + noise
signal1 = sin(sin_freq*2*pi*t)+1/10*randn(size(t));

% signal 2 is in phase
signal2 = sin(sin_freq*2*pi*t)+1/10*randn(size(t));

% signal 3 is 180 out of phase
signal3 = sin(pi+sin_freq*2*pi*t)+1/10*randn(size(t));

figure('Position',[0 0 800 300],'Color',[1 1 1])

bar_width = .3;

% signal 1 and 2
subplot(4,4,1), hold on
plot(t,signal1,'k')
title('V_1')
subplot(4,8,4)
bar(1,mean(signal1.^2),bar_width,'w')
title('mn(V_1^2)')

subplot(4,4,5), hold on
plot(t,signal2,'k')
title('V_2')
subplot(4,8,12)
bar(1,mean(signal2.^2),bar_width,'w')
title('mn(V_2^2)')

subplot(4,4,9), hold on
plot(t,signal1 + signal2,'k')
title('V_1 + V_2')
subplot(4,8,19)
bar(1,mean((signal1 + signal2).^2),bar_width,'w')
title('mn((V_1 + V_2)^2)')

subplot(4,8,20)
bar(1,mean(signal1.^2) + mean(signal2.^2),bar_width,'w')
title('mn(V_1^2) + mn(V_2^2)')


% signal 1 and 3
subplot(4,4,3), hold on
plot(t,signal1,'k')
title('V_1')
subplot(4,8,8)
bar(1,mean(signal1.^2),bar_width,'w')
title('mn(V_1^2)')

subplot(4,4,7), hold on
plot(t,signal3,'k')
title('V_3')
subplot(4,8,16)
bar(1,mean(signal3.^2),bar_width,'w')
title('mn(V_3^2)')

subplot(4,4,11), hold on
plot(t,signal1 + signal3,'k')
title('V_1 + V_3')
subplot(4,8,23)
bar(1,mean((signal1 + signal3).^2),bar_width,'w')
title('mn((V_1 + V_3)^2)')

subplot(4,8,24)
bar(1,mean(signal1.^2) + mean(signal3.^2),bar_width,'w')
title('mn(V_1^2) + mn(V_3^2)')

for k = [1:2:16]
    subplot(4,4,k)
    ylim([-3 3])
end


for k = [4 12 19 20 8 16 23 24]
    subplot(4,8,k)
    ylim([-0.1 2.1])
    set(gca,'XTick',[],'YTick',[0 1 2])
    box off
end

set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/ns_ConceptMethodsFig01'])
% print('-depsc','-r300',['./figures/ns_ConceptMethodsFig01'])


%%
% This script generates pannels for Fig 2 from Hermes et al:
%
% This script generates a set of signals with varying correlation
% (coherence), variance, and mean. It shows how the sum of the power, the
% power of the sum and the cross-power vary as a function of these three
% factors.
%
% JW and DH 2017

%% start with one example:

m = 2;      % mean
v = 2;      % variance
s = .1;     % coherence
n = 5;      % number of neurons
T = 1000;   % number of time points
dt = 1/T;

% covariance matrix
sigma = v * (eye(n) + (1-eye(n)) * s);

% data
x = mvnrnd(ones(T,n) * m, sigma);

% observed total cross power
% this should have a time factor added for integration period
total_cp  = x' * x * dt; 

% partitioned into diagonal (power) and off-diagonal ('cross power')
pw  = sum(diag(total_cp));
cpw = sum(sum(total_cp .* (1-eye(n))));

disp([pw cpw])

% predicted
pw_pred  = (m^2 + v)*n;
cpw_pred = (m^2 + v*s)*(n^2-n);
disp([pw_pred cpw_pred])

disp('note that these are similar only if data are 1 sec long or add the time period to total_cp')

%% Make the figure:
%% 
%% Plot total power (LFP), power (BOLD) and cross-power
%%

figure('Position',[0 0 700 250])

fs = 10;    % font size
n = 10;     % number of neurons
[v, s] = meshgrid(linspace(0,1,10));

for ii = 1:2
    
    if ii == 1, m = 0; else m = 1; end
    
    pw_pred  = (m^2 + v)*n;
    cpw_pred = (m^2 + v.*s)*(n^2-n);
    
    subplot(1,3,1),set(gca, 'FontSize', fs);
    total_pw = cpw_pred+pw_pred;
    if ii==2
        mesh(v, s, total_pw, total_pw - min(total_pw(:))),  hold on
    else
        surf(v, s, total_pw, total_pw - min(total_pw(:))),  hold on
    end
    title('Power of sum (LFP)')
    zlim([0 200])
    
    subplot(1,3,2), set(gca, 'FontSize', fs);
    if ii==2
        sH = mesh(v, s, pw_pred, pw_pred - min(pw_pred(:)));  hold on
    else
        sH = surf(v, s, pw_pred, pw_pred - min(pw_pred(:)));  hold on
    end
    set(sH, 'FaceLighting', 'gouraud')
    title('Sum of power (BOLD)')
    zlim([0 200]) 
    
    subplot(1,3,3),set(gca, 'FontSize', fs);
    if ii==2
        mesh(v, s, cpw_pred, cpw_pred-min(cpw_pred(:))), hold on
    else
        surf(v, s, cpw_pred, cpw_pred-min(cpw_pred(:))), hold on
    end
    title('Sum of cross Power')
    zlim([0 200])
end

for k = 1:3
    subplot(1,3,k)
    h=xlabel('Variance (\sigma^2)'); set(h, 'rotation', 10)
    h=ylabel('Correlation (\rho)');set(h, 'rotation', -30)
    alpha(.8)
    view(-30,10)
    set(gca,'ZTick',[0:20:200])
end

colormap parula

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['../figures/analytic/2D_LFP_BOLD_CPW'])
% print('-dpng','-r300',['../figures/analytic/2D_LFP_BOLD_CPW'])


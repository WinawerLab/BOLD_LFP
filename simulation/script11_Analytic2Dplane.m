m = 2;      % mean
v = 2;      % variance
s = .1;     % coherence
n = 5;      % number of neurons
T = 10000;  % number of time points
dt = 1/T;

% covariance matrix
sigma = v * (eye(n) + (1-eye(n)) * s);

% data
x = mvnrnd(ones(T,n) * m, sigma);

% observed total cross power
total_cp  = x' * x * dt;

% partitioned into diagonal (power) and off-diagonal ('cross power')
pw  = sum(diag(total_cp));
cpw = sum(sum(total_cp .* (1-eye(n))));

disp([pw cpw])

% predicted
pw_pred  = (m^2 + v)*n;
cpw_pred = (m^2 + v*s)*(n^2-n);
disp([pw_pred cpw_pred])

%%
%% PLOT NICE 2D function
%%

fH = figure; clf; pos=get(fH, 'Position');
pos([3 4]) = [1200 400]; set(fH, 'Position', pos);

fs = 18;    % font size
n = 10;     % number of neurons
[v, s] = meshgrid(linspace(0,1, 20));

for ii = 1:2
    
    if ii == 1, m = 0; else m = 1; end
    
    pw_pred  = (m^2 + v)*n;
    cpw_pred = (m^2 + v.*s)*(n^2-n);
    
    figure(1);
    subplot(1,3,1), set(gca, 'FontSize', fs);
    sH = surf(v, s, pw_pred, pw_pred - min(pw_pred(:)));  hold on
    set(sH, 'FaceLighting', 'gouraud')
    h=xlabel('Variance (\sigma^2)'); set(h, 'rotation', 30)
    h=ylabel('Correlation (\rho)');set(h, 'rotation', -45)
    title('Power (BOLD)')
    
    subplot(1,3,2),set(gca, 'FontSize', fs);
    surf(v, s, cpw_pred, cpw_pred-min(cpw_pred(:))), hold on
    h=xlabel('Variance (\sigma^2)'); set(h, 'rotation', 30)
    h=ylabel('Correlation (\rho)');set(h, 'rotation', -45)
    title('Cross Power')
    
    subplot(1,3,3),set(gca, 'FontSize', fs);
    total_pw = cpw_pred+pw_pred;
    surf(v, s, total_pw, total_pw - min(total_pw(:))),  hold on
    h=xlabel('Variance (\sigma^2)'); set(h, 'rotation', 30)
    h=ylabel('Correlation (\rho)');set(h, 'rotation', -45)
    title('Total Power (LFP)')
    
    %     subplot(2,3,4), set(gca, 'ColorOrder', jet(size(v,1)), 'FontSize', fs); hold on;
    %     plot(v(1,:), pw_pred')
    %     xlabel('Variance');
    %     title('Power (BOLD)')
    %
    %     subplot(2,3,5), set(gca, 'ColorOrder', jet(size(v,1)), 'FontSize', fs); hold on;
    %     plot(v(1,:),cpw_pred')
    %     xlabel('Variance');
    %     title('Cross Power')
    %
    %     subplot(2,3,6), set(gca, 'ColorOrder', jet(size(v,1)), 'FontSize', fs); hold on;
    %     plot(v(1,:), cpw_pred'+pw_pred')
    %     xlabel('Variance');
    %     title('Total Power (LFP)')
    
end
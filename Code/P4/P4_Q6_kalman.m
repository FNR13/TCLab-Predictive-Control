% P4 Q6 - augmented Kalman filter, tested in OPEN LOOP (no MPC).
%
% The simulator has the perturbed c1 (same +12.5 % input-bias mismatch as
% Q4): the filter has a non-zero disturbance to estimate.
% The applied control is the constant feedforward u = u_ss (i.e. Du = 0):
% no feedback loop, so we observe the filter behaviour in isolation.
%
% Three values of the disturbance-model variance delta_E are compared:
%   - small  (1e-6): filter trusts the constant-disturbance assumption,
%                    estimates d very slowly;
%   - medium (1e-3): the chosen default for Q7;
%   - large  (1e+0): the filter reacts fast but is noisy.
%
% Initial state estimate has a built-in error equivalent to an output
% offset of about +5 deg C, to make the convergence visible.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
common_dir = fullfile(script_dir, '..', 'Common');
addpath(common_dir);
out_dir    = fullfile(script_dir, '..', '..', '..', 'Imagens');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

set(0,'DefaultFigureColor','w','DefaultAxesColor','w', ...
      'DefaultAxesXColor','k','DefaultAxesYColor','k', ...
      'DefaultAxesGridColor',[0.15 0.15 0.15],'DefaultTextColor','k', ...
      'DefaultAxesFontSize',11,'DefaultLegendColor','w', ...
      'DefaultLegendTextColor','k','DefaultLegendEdgeColor','k', ...
      'DefaultFigureInvertHardcopy','off');
save_white = @(fig, name) exportgraphics(fig, fullfile(out_dir,name), ...
                                         'BackgroundColor','white','Resolution',150);

% ----- model -------------------------------------------------------------
load(fullfile(common_dir,'singleheater_model.mat'), ...
     'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n     = size(A,1);
e_std = sqrt(e_var);

x_ss   = [eye(n)-A; C] \ [B*u_ss; y_ss];
c1_nom = (eye(n)-A)*x_ss - B*u_ss;
c2     = y_ss - C*x_ss;

% Same mismatch chosen in Q4
d_real  = 0.50 * u_ss;
c1_pert = c1_nom + B * d_real;

% Plant in this Q6 uses the PERTURBED c1 (so there is a real d to estimate)
h1 = @(x,u) A*x + B*u + Ke*e_std*randn + c1_pert;
T1C = @(x) C*x + e_std*randn + c2;

% ----- simulation parameters --------------------------------------------
T = 1500; N = T/Ts; t = (0:N-1)*Ts;

% Start at ambient temperature
Dx0  = (eye(n)-A)\(B*(0 - u_ss));
x0   = Dx0 + x_ss;

% Initial state ESTIMATE: add an error so Dx_hat(0) corresponds to a
% measurement about +5 C above the real one. C * Dx_err = +5 means we
% need Dx_err such that C * Dx_err = 5. Pick the minimum-norm one.
Dx_err = pinv(C) * 5;       % gives C*Dx_err = 5

% Three filter variances to sweep
delta_E_list = [1e-6, 1e-3, 1e+0];
labels       = {'\delta_E = 10^{-6}', '\delta_E = 10^{-3}', '\delta_E = 10^{+0}'};
colors       = lines(numel(delta_E_list));

% Storage
Y     = nan(N, numel(delta_E_list));
Yhat  = nan(N, numel(delta_E_list));
Dhat  = nan(N, numel(delta_E_list));

fprintf('--- P4 Q6: augmented Kalman, open loop, u = u_ss ---\n');
fprintf('real disturbance d_real = %.4f (= %.1f %% of u_ss)\n', d_real, 100*d_real/u_ss);
fprintf('\n%-22s %-14s %-14s %-14s\n', 'delta_E', 'd_hat (final)', 'y err (final)', 't_50%% [s]');

for jj = 1:numel(delta_E_list)
    delta_E = delta_E_list(jj);
    [A_d, B_d, C_d, L] = kalman_augmented(A, B, C, Ke, e_var, delta_E);

    rng(0,'twister');
    x = nan(n,N+1); x(:,1) = x0;
    xd_hat = [Dx0 + Dx_err; 0];      % initial estimate with +5 C error
    Du_prev = 0;                     % we apply u = u_ss, so Du = 0

    for k = 1:N
        % measurement
        y_k  = T1C(x(:,k));
        Dy_k = y_k - y_ss;
        Y(k,jj) = y_k;

        % Kalman: predict + correct
        xd_minus = A_d * xd_hat + B_d * Du_prev;
        xd_hat   = xd_minus + L * (Dy_k - C_d * xd_minus);

        Yhat(k,jj) = C * xd_hat(1:n) + y_ss;   % output as estimated
        Dhat(k,jj) = xd_hat(n+1);

        % open-loop control:  u = u_ss
        u_k = u_ss;
        x(:,k+1) = h1(x(:,k), u_k);

        Du_prev = 0;
    end

    % metrics
    d_final  = Dhat(end,jj);
    err_final = Yhat(end,jj) - Y(end,jj);
    % 50 % settling: first k where |d_hat(k) - d_real| < 0.5 * d_real
    settled = find(abs(Dhat(:,jj) - d_real) < 0.5*d_real, 1, 'first');
    if isempty(settled), t_settle = NaN; else, t_settle = t(settled); end

    fprintf('%-22s %-14.4f %-14.4f %-14.1f\n', ...
            labels{jj}, d_final, err_final, t_settle);
end

% ----- figures ----------------------------------------------------------
% (a) measured y and estimated y (single delta_E = 1e-3)
jj_pick = 2;
fig = figure('Position',[100 100 950 420],'Color','w');
hold on; grid on; box on;
plot(t, Y(:,jj_pick), 'k', 'LineWidth',1.1, 'DisplayName','y (measured)');
plot(t, Yhat(:,jj_pick), 'r--', 'LineWidth',1.2, 'DisplayName','y\^ (estimated)');
yline(y_ss,'k:','y_{ss}');
title(sprintf('Q6 - Kalman filter open loop (u = u_{ss}), %s', labels{jj_pick}));
xlabel('t [s]'); ylabel('y [°C]'); legend('Location','best');
save_white(fig,'Q6_y_vs_yhat.png'); close(fig);

% (b) d_hat for the three values of delta_E
fig = figure('Position',[100 100 950 420],'Color','w');
hold on; grid on; box on;
for jj = 1:numel(delta_E_list)
    plot(t, Dhat(:,jj), 'LineWidth',1.2, 'Color',colors(jj,:), ...
         'DisplayName',labels{jj});
end
yline(d_real,'k--',sprintf('d_{real} = %.3f', d_real));
title('Q6 - Estimated input disturbance d\^ for different \delta_E');
xlabel('t [s]'); ylabel('d\^');
legend('Location','best');
save_white(fig,'Q6_dhat_vs_deltaE.png'); close(fig);

% (c) estimation error y_hat - y
fig = figure('Position',[100 100 950 420],'Color','w');
hold on; grid on; box on;
for jj = 1:numel(delta_E_list)
    plot(t, Yhat(:,jj) - Y(:,jj), 'LineWidth',1.0, 'Color',colors(jj,:), ...
         'DisplayName',labels{jj});
end
yline(0,'k:');
title('Q6 - Output estimation error  y\^ - y');
xlabel('t [s]'); ylabel('y\^ - y [°C]');
legend('Location','best');
save_white(fig,'Q6_yerr_vs_deltaE.png'); close(fig);

fprintf('Figures written to %s\n', out_dir);

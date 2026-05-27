% P4 Q7 - full closed loop: MPC + Kalman + tracking + soft output cap.
%
% Reference profile (absolute): 50 -> 40 -> 60 -> 45 deg C.
%   60 C is intentionally above the safety cap (55 C) so the soft cap
%   must clip it.
% The simulator runs with the perturbed c1 (input bias d_real = 12.5 %),
% which the augmented Kalman filter must estimate to remove the offset
% that was visible in Q4 (no integral action there).
%
% Compares two configurations:
%   (a) WITH disturbance estimate fed into the feedforward target
%       (proper zero-offset MPC, the controller designed for the project);
%   (b) WITHOUT disturbance estimate (d_hat clamped to 0 in
%       mpc_steady_state) -- to illustrate that the Kalman estimate is
%       what produces the zero offset.

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

% Same +12.5 % input bias as Q4 (kept on, per the project guide)
d_real  = 0.50 * u_ss;
c1_pert = c1_nom + B * d_real;

h1   = @(x,u) A*x + B*u + Ke*e_std*randn + c1_pert;
T1C  = @(x)   C*x + e_std*randn + c2;

% ----- design ------------------------------------------------------------
H      = 20;
R      = 1e-2;
y_max  = 55;
alpha  = 1e4;
deltaE = 1e-1;        % chosen from Q6: faster than 1e-3, less noisy than 1e0

[A_d, B_d, C_d, L] = kalman_augmented(A, B, C, Ke, e_var, deltaE);

% ----- reference profile -------------------------------------------------
r_levels = [50, 40, 60, 45];          % C
T_step   = 600;                       % s per level
N_step   = T_step / Ts;
N        = N_step * numel(r_levels);
T        = N * Ts;
t        = (0:N-1)*Ts;
r        = zeros(1, N);
for jj = 1:numel(r_levels)
    r((jj-1)*N_step + (1:N_step)) = r_levels(jj);
end

% ----- initial condition (ambient) --------------------------------------
Dx0  = (eye(n)-A)\(B*(0 - u_ss));
x0   = Dx0 + x_ss;

% Run both configurations
configs = {'with d_hat', 'without d_hat'};
use_d   = [true, false];

results = struct();
for cfg = 1:2
    rng(0,'twister');
    x  = nan(n, N+1); x(:,1) = x0;
    y  = nan(1,N);    u   = nan(1,N);
    Dxh   = nan(n,N); dh  = nan(1,N);
    exitf = nan(1,N); tq  = nan(1,N);

    xd_hat   = [zeros(n,1); 0];    % start with zero estimate
    Du_prev  = 0;

    for k = 1:N
        % measurement
        y(k)  = T1C(x(:,k));
        Dy_k  = y(k) - y_ss;

        % Kalman: predict + correct
        xd_minus = A_d * xd_hat + B_d * Du_prev;
        xd_hat   = xd_minus + L * (Dy_k - C_d * xd_minus);
        Dx_hat   = xd_hat(1:n);
        d_hat    = xd_hat(n+1);
        Dxh(:,k) = Dx_hat;
        dh(k)    = d_hat;

        % tracking target
        Dr = r(k) - y_ss;
        if use_d(cfg)
            [Dx_bar, Du_bar] = mpc_steady_state(A,B,C, Dr, d_hat);
        else
            [Dx_bar, Du_bar] = mpc_steady_state(A,B,C, Dr, 0);
        end

        % MPC options (shifted bounds and shifted soft cap)
        opts = struct('formulation','dense', ...
                      'u_min',  0  - u_ss - Du_bar, ...
                      'u_max', 100 - u_ss - Du_bar, ...
                      'y_max', (y_max - y_ss) - Dr, ...
                      'soft',  true, 'alpha', alpha);

        delta_x = Dx_hat - Dx_bar;

        tic;
        [delta_u, info] = mpc_solve(delta_x, A, B, C, H, R, opts);
        tq(k) = toc;
        exitf(k) = info.exitflag;

        Du_k   = delta_u + Du_bar;
        u(k)   = u_ss + Du_k;
        x(:,k+1) = h1(x(:,k), u(k));
        Du_prev = Du_k;
    end

    results(cfg).y    = y;    results(cfg).u    = u;
    results(cfg).dh   = dh;   results(cfg).ef   = exitf;
    results(cfg).tq   = tq;
end

% ----- offset metrics per reference segment -----------------------------
fprintf('--- P4 Q7: full closed loop, reference profile [50, 40, 60, 45] C ---\n');
fprintf('Augmented Kalman: delta_E = %g, ||L|| = %.3f\n', deltaE, norm(L));
fprintf('Real disturbance d_real = %.2f, soft cap y_max = %.1f C\n\n', d_real, y_max);

for cfg = 1:2
    fprintf('  -- %s --\n', configs{cfg});
    fprintf('    %-6s  %-9s  %-9s  %-9s\n', 'level', 'r [C]', 'y_mean','offset');
    for jj = 1:numel(r_levels)
        idx = (jj-1)*N_step + (1:N_step);
        % take only the last 40 % of the segment (after settling)
        tail = idx(round(0.6*N_step):end);
        y_mean = mean(results(cfg).y(tail));
        target = min(r_levels(jj), y_max);   % effective reference under cap
        fprintf('    %-6d  %-9.2f  %-9.2f  %+8.3f\n', ...
                jj, r_levels(jj), y_mean, y_mean - target);
    end
    fprintf('    final d_hat = %+8.3f (real %.3f)\n', results(cfg).dh(end), d_real);
    fprintf('    mean tQP = %.2f ms,  exitflag != 1 count = %d\n\n', ...
            mean(results(cfg).tq)*1e3, sum(results(cfg).ef ~= 1));
end

% ----- figure (a): comparison output / control --------------------------
fig = figure('Position',[100 100 1000 750],'Color','w');

subplot(3,1,1); hold on; grid on; box on;
title('Q7 - closed loop: reference tracking with safety cap (perturbed plant)');
plot(t, results(1).y, 'b', 'LineWidth',1.2, 'DisplayName','y, with d\^');
plot(t, results(2).y, 'r', 'LineWidth',1.2, 'DisplayName','y, without d\^');
stairs(t, r, 'k--', 'LineWidth',1.0, 'DisplayName','reference r');
yline(y_max,'r-','y_{max}');
ylabel('y [°C]'); legend('Location','best');

subplot(3,1,2); hold on; grid on; box on;
stairs(t, results(1).u, 'b', 'LineWidth',1.0, 'DisplayName','u, with d\^');
stairs(t, results(2).u, 'r', 'LineWidth',1.0, 'DisplayName','u, without d\^');
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); legend('Location','best');

subplot(3,1,3); hold on; grid on; box on;
plot(t, results(1).dh, 'b', 'LineWidth',1.0, 'DisplayName','d\^, with d\^');
plot(t, results(2).dh, 'r', 'LineWidth',1.0, 'DisplayName','d\^, without d\^');
yline(d_real,'k--', sprintf('d_{real} = %.2f', d_real));
ylim([-20, 20]);    % crop the transient burst at start-up
ylabel('d\^'); xlabel('t [s]'); legend('Location','best');

save_white(fig,'Q7_full_closed_loop.png'); close(fig);
fprintf('Figure written: %s\n', fullfile(out_dir,'Q7_full_closed_loop.png'));

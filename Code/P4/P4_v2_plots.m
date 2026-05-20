% P4 v2 plots - refined figures for the final report.
%
% Re-runs the same simulations as P4_Q2..Q7 (same rng(0) seed -> identical
% traces and identical numbers) and writes a "_v2" set of PNGs to
% ECPD/Imagens, WITHOUT touching the original files.
%
% Changes vs the original figures:
%   * Q2 (sweepR, sweepH): zoom the time axis to 0-500 s with xlim. The
%     simulation still runs 1500 s so the metrics behind Table 2 are
%     unchanged; only the displayed window is shorter. The heatmap is NOT
%     regenerated.
%   * Q3, Q4: same 0-500 s zoom.
%   * Q6 (dhat, yerr): zoom to 0-750 s. The y-vs-yhat figure is left as is
%     (no v2).
%   * Q3, Q4, Q5, Q7: reference / limit lines recoloured. Input limit
%     lines u = 0 % and u = 100 % are hidden from the legend; interpretive
%     lines (reference, steady state, safety cap, true disturbance) are
%     given explicit English labels and distinct colours.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
common_dir = fullfile(script_dir, '..', 'Common');
addpath(common_dir);
out_dir    = fullfile(script_dir, '..', '..', '..', 'Imagens');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

% white-background defaults
set(0,'DefaultFigureColor','w','DefaultAxesColor','w', ...
      'DefaultAxesXColor','k','DefaultAxesYColor','k', ...
      'DefaultAxesGridColor',[0.15 0.15 0.15],'DefaultTextColor','k', ...
      'DefaultAxesFontSize',11,'DefaultLegendColor','w', ...
      'DefaultLegendTextColor','k','DefaultLegendEdgeColor','k', ...
      'DefaultColorbarColor','k','DefaultFigureInvertHardcopy','off');
save_v2 = @(fig, name) exportgraphics(fig, fullfile(out_dir,name), ...
                                      'BackgroundColor','white','Resolution',150);

% neutral colours for auxiliary / interpretive lines
COL_LIMIT = [0.35 0.35 0.35];   % grey  - input bounds (hidden from legend)
COL_REF   = [0.00 0.00 0.00];   % black - reference / steady state
COL_CAP   = [0.45 0.10 0.70];   % purple - safety cap
COL_DIST  = [0.00 0.60 0.30];   % green  - true disturbance

% ----- shared model ------------------------------------------------------
load(fullfile(common_dir,'singleheater_model.mat'), ...
     'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n     = size(A,1);
e_std = sqrt(e_var);

x_ss   = [eye(n)-A; C] \ [B*u_ss; y_ss];
c1_nom = (eye(n)-A)*x_ss - B*u_ss;
c2     = y_ss - C*x_ss;

h1_make = @(c1) @(x,u) A*x + B*u + Ke*e_std*randn + c1;
T1C     = @(x)   C*x + e_std*randn + c2;

H = 20; R = 1e-2;
T = 1500;  N = T/Ts;
Dx0_amb = (eye(n)-A)\(B*(0 - u_ss));
x0_amb  = Dx0_amb + x_ss;

% =========================================================================
% Q2 - sweep H and R, unconstrained  (display window 0-500 s)
% =========================================================================
fprintf('Q2 v2...\n');

H_grid = [5 10 15 20 25 30];
R_grid = [1e-3, 1e-2, 1e-1, 1, 10];
nH = numel(H_grid); nR = numel(R_grid);
trace_y = cell(nH,nR); trace_u = cell(nH,nR); trace_t = cell(nH,nR);

for ih = 1:nH
    for ir = 1:nR
        H2 = H_grid(ih); R2 = R_grid(ir);
        rng(0,'twister');
        x = nan(n,N+1); x(:,1) = x0_amb;
        y = nan(1,N); u = nan(1,N); t = (0:N-1)*Ts;
        h1 = h1_make(c1_nom);
        opts = struct('formulation','dense');
        for k = 1:N
            y(k) = T1C(x(:,k));
            Dx_k = x(:,k) - x_ss;
            Du_k = mpc_solve(Dx_k, A, B, C, H2, R2, opts);
            u(k) = u_ss + Du_k;
            x(:,k+1) = h1(x(:,k), u(k));
        end
        trace_y{ih,ir} = y; trace_u{ih,ir} = u; trace_t{ih,ir} = t;
    end
end

ih_plot = find(H_grid == 20);
ir_plot = find(abs(R_grid - 1e-2) < 1e-9);
XLIM_Q2 = [0 500];

% Q2 sweep R at H=20
fig = figure('Position',[100 100 900 600],'Color','w');
subplot(2,1,1); hold on; grid on; box on;
title('Q2 sweep of R at H = 20 (unconstrained)');
for ir = 1:nR
    plot(trace_t{ih_plot,ir}, trace_y{ih_plot,ir}, ...
         'LineWidth',1.0, 'DisplayName', sprintf('R = %g', R_grid(ir)));
end
yline(y_ss,'--','Color',COL_REF,'Label','Steady state', ...
      'LineWidth',1.0,'HandleVisibility','off');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
xlim(XLIM_Q2);
subplot(2,1,2); hold on; grid on; box on;
for ir = 1:nR
    stairs(trace_t{ih_plot,ir}, trace_u{ih_plot,ir}, 'LineWidth',1.0);
end
yline(0,'--','Color',COL_LIMIT,'HandleVisibility','off');
yline(100,'--','Color',COL_LIMIT,'HandleVisibility','off');
ylabel('u [%]'); xlabel('t [s]');
xlim(XLIM_Q2);
save_v2(fig, 'Q2_sweepR_H20_v2.png'); close(fig);

% Q2 sweep H at R=0.01
fig = figure('Position',[100 100 900 600],'Color','w');
subplot(2,1,1); hold on; grid on; box on;
title('Q2 sweep of H at R = 0.01 (unconstrained)');
for ih = 1:nH
    plot(trace_t{ih,ir_plot}, trace_y{ih,ir_plot}, ...
         'LineWidth',1.0, 'DisplayName', sprintf('H = %d', H_grid(ih)));
end
yline(y_ss,'--','Color',COL_REF,'Label','Steady state', ...
      'LineWidth',1.0,'HandleVisibility','off');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
xlim(XLIM_Q2);
subplot(2,1,2); hold on; grid on; box on;
for ih = 1:nH
    stairs(trace_t{ih,ir_plot}, trace_u{ih,ir_plot}, 'LineWidth',1.0);
end
yline(0,'--','Color',COL_LIMIT,'HandleVisibility','off');
yline(100,'--','Color',COL_LIMIT,'HandleVisibility','off');
ylabel('u [%]'); xlabel('t [s]');
xlim(XLIM_Q2);
save_v2(fig, 'Q2_sweepH_R0.01_v2.png'); close(fig);

% =========================================================================
% Q3 - bounded vs unbounded  (display window 0-500 s)
% =========================================================================
fprintf('Q3 v2...\n');

Du_min = 0   - u_ss;
Du_max = 100 - u_ss;
opts_unb = struct('formulation','dense');
opts_box = struct('formulation','dense','u_min',Du_min,'u_max',Du_max);

[t_u,y_u,u_u] = run_loop_reg(A,B,C, h1_make(c1_nom), T1C, x0_amb, x_ss, u_ss, ...
                             N, Ts, H, R, opts_unb);
[t_b,y_b,u_b] = run_loop_reg(A,B,C, h1_make(c1_nom), T1C, x0_amb, x_ss, u_ss, ...
                             N, Ts, H, R, opts_box);

fig = figure('Position',[100 100 950 650],'Color','w');
subplot(2,1,1); hold on; grid on; box on;
title('Q3 - closed loop, H = 20, R = 0.01');
plot(t_u, y_u, 'LineWidth',1.2,'DisplayName','y, unbounded');
plot(t_b, y_b, 'LineWidth',1.2,'DisplayName','y, bounded [0,100]');
yline(y_ss,'--','Color',COL_REF,'Label','Steady state', ...
      'LineWidth',1.0,'HandleVisibility','off');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
xlim([0 500]);
subplot(2,1,2); hold on; grid on; box on;
stairs(t_u, u_u, 'LineWidth',1.0,'DisplayName','u, unbounded');
stairs(t_b, u_b, 'LineWidth',1.2,'DisplayName','u, bounded');
yline(0,'--','Color',COL_LIMIT,'HandleVisibility','off');
yline(100,'--','Color',COL_LIMIT,'HandleVisibility','off');
ylabel('u [%]'); xlabel('t [s]'); legend('Location','best');
xlim([0 500]);
save_v2(fig, 'Q3_bounded_vs_unbounded_v2.png'); close(fig);

% =========================================================================
% Q4 - tracking nominal vs perturbed  (display window 0-500 s)
% =========================================================================
fprintf('Q4 v2...\n');

Dr     = 5;  r_abs = y_ss + Dr;
[Dx_bar, Du_bar] = mpc_steady_state(A, B, C, Dr, 0);
opts_t = struct('formulation','dense', ...
                'u_min', Du_min - Du_bar, ...
                'u_max', Du_max - Du_bar);

d_real  = 0.50 * u_ss;
c1_pert = c1_nom + B * d_real;

[t_a,y_a,u_a] = run_loop_track(A,B,C, h1_make(c1_nom), T1C, x0_amb, ...
                               x_ss, u_ss, Dx_bar, Du_bar, N, Ts, H, R, opts_t);
[t_b4,y_b4,u_b4] = run_loop_track(A,B,C, h1_make(c1_pert), T1C, x0_amb, ...
                                  x_ss, u_ss, Dx_bar, Du_bar, N, Ts, H, R, opts_t);

fig = figure('Position',[100 100 950 700],'Color','w');
subplot(2,1,1); hold on; grid on; box on;
title(sprintf('Q4 - tracking r = %.2f C, H = 20, R = 0.01', r_abs));
plot(t_a, y_a, 'LineWidth',1.2,'DisplayName','y, nominal');
plot(t_b4, y_b4, 'LineWidth',1.2,'DisplayName','y, perturbed');
yline(r_abs,'--','Color',COL_REF,'Label','Reference', ...
      'LineWidth',1.1,'HandleVisibility','off');
yline(y_ss, ':','Color',COL_CAP,'Label','Steady state', ...
      'LineWidth',1.1,'HandleVisibility','off');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
xlim([0 500]);
subplot(2,1,2); hold on; grid on; box on;
stairs(t_a, u_a, 'LineWidth',1.0,'DisplayName','u, nominal');
stairs(t_b4, u_b4, 'LineWidth',1.0,'DisplayName','u, perturbed');
yline(0,'--','Color',COL_LIMIT,'HandleVisibility','off');
yline(100,'--','Color',COL_LIMIT,'HandleVisibility','off');
ylabel('u [%]'); xlabel('t [s]'); legend('Location','best');
xlim([0 500]);
save_v2(fig, 'Q4_tracking_nominal_vs_perturbed_v2.png'); close(fig);

% =========================================================================
% Q5 - hard vs soft output cap  (recoloured lines, full time axis)
% =========================================================================
fprintf('Q5 v2...\n');

r_abs5 = 60; Dr5 = r_abs5 - y_ss;
[Dx_bar5, Du_bar5] = mpc_steady_state(A, B, C, Dr5, 0);
y_max  = 55;
dy_max = (y_max - y_ss) - Dr5;
opts_hard = struct('formulation','dense', ...
                   'u_min', 0 - u_ss - Du_bar5, ...
                   'u_max', 100 - u_ss - Du_bar5, ...
                   'y_max', dy_max, 'soft', false);
opts_soft = opts_hard; opts_soft.soft = true; opts_soft.alpha = 1e4;

u_hot = 70;
Dx0_hot = (eye(n)-A)\(B*(u_hot - u_ss));
x0_hot  = Dx0_hot + x_ss;

[t_h,y_h,u_h] = run_loop_track(A,B,C, h1_make(c1_nom), T1C, x0_hot, ...
                               x_ss, u_ss, Dx_bar5, Du_bar5, N, Ts, H, R, opts_hard);
[t_s,y_s,u_s] = run_loop_track(A,B,C, h1_make(c1_nom), T1C, x0_hot, ...
                               x_ss, u_ss, Dx_bar5, Du_bar5, N, Ts, H, R, opts_soft);

fig = figure('Position',[100 100 950 700],'Color','w');
subplot(2,1,1); hold on; grid on; box on;
title(sprintf('Q5 - reference r = %.1f C with cap y_{max} = %.1f C', r_abs5, y_max));
plot(t_h, y_h, 'LineWidth',1.2,'DisplayName','y, hard cap');
plot(t_s, y_s, 'LineWidth',1.2,'DisplayName','y, soft cap');
yline(r_abs5,'--','Color',COL_REF,'Label','Reference', ...
      'LineWidth',1.1,'HandleVisibility','off');
yline(y_max, '-','Color',COL_CAP,'Label','Safety cap', ...
      'LineWidth',1.3,'HandleVisibility','off');
yline(y_ss,  ':','Color',COL_LIMIT,'Label','Steady state', ...
      'LineWidth',1.0,'HandleVisibility','off');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
subplot(2,1,2); hold on; grid on; box on;
stairs(t_h, u_h, 'LineWidth',1.0,'DisplayName','u, hard cap');
stairs(t_s, u_s, 'LineWidth',1.0,'DisplayName','u, soft cap');
yline(0,'--','Color',COL_LIMIT,'HandleVisibility','off');
yline(100,'--','Color',COL_LIMIT,'HandleVisibility','off');
ylabel('u [%]'); xlabel('t [s]'); legend('Location','best');
save_v2(fig, 'Q5_hard_vs_soft_v2.png'); close(fig);

% =========================================================================
% Q6 - augmented Kalman, dhat and yerr  (display window 0-750 s)
% =========================================================================
fprintf('Q6 v2...\n');

d_real6 = 0.50 * u_ss;
c1_p6   = c1_nom + B * d_real6;
h1_6    = h1_make(c1_p6);

Dx_err = pinv(C) * 5;
delta_E_list = [1e-6, 1e-3, 1e+0];
labels6      = {'\delta_E = 10^{-6}', '\delta_E = 10^{-3}', '\delta_E = 10^{+0}'};
colors6      = lines(numel(delta_E_list));

t6 = (0:N-1)*Ts;
Y6    = nan(N, numel(delta_E_list));
Yhat6 = nan(N, numel(delta_E_list));
Dhat6 = nan(N, numel(delta_E_list));

for jj = 1:numel(delta_E_list)
    [A_d, B_d, C_d, L] = kalman_augmented(A, B, C, Ke, e_var, delta_E_list(jj));
    rng(0,'twister');
    x = nan(n,N+1); x(:,1) = x0_amb;
    xd_hat = [Dx0_amb + Dx_err; 0];
    Du_prev = 0;
    for k = 1:N
        y_k  = T1C(x(:,k));
        Dy_k = y_k - y_ss;
        Y6(k,jj) = y_k;
        xd_minus = A_d * xd_hat + B_d * Du_prev;
        xd_hat   = xd_minus + L * (Dy_k - C_d * xd_minus);
        Yhat6(k,jj) = C * xd_hat(1:n) + y_ss;
        Dhat6(k,jj) = xd_hat(n+1);
        x(:,k+1) = h1_6(x(:,k), u_ss);
        Du_prev  = 0;
    end
end

XLIM_Q6 = [0 750];

% d_hat for the three delta_E
fig = figure('Position',[100 100 950 420],'Color','w');
hold on; grid on; box on;
for jj = 1:numel(delta_E_list)
    plot(t6, Dhat6(:,jj), 'LineWidth',1.2, 'Color',colors6(jj,:), ...
         'DisplayName',labels6{jj});
end
yline(d_real6,'--','Color',COL_DIST, ...
      'Label',sprintf('True disturbance (%.2f)', d_real6), ...
      'LineWidth',1.3,'HandleVisibility','off');
title('Q6 - Estimated input disturbance d\^ for different \delta_E');
xlabel('t [s]'); ylabel('d\^'); legend('Location','best');
xlim(XLIM_Q6);
save_v2(fig,'Q6_dhat_vs_deltaE_v2.png'); close(fig);

% estimation error y_hat - y
fig = figure('Position',[100 100 950 420],'Color','w');
hold on; grid on; box on;
for jj = 1:numel(delta_E_list)
    plot(t6, Yhat6(:,jj) - Y6(:,jj), 'LineWidth',1.0, 'Color',colors6(jj,:), ...
         'DisplayName',labels6{jj});
end
yline(0,':','Color',COL_LIMIT,'LineWidth',1.0,'HandleVisibility','off');
title('Q6 - Output estimation error  y\^ - y');
xlabel('t [s]'); ylabel('y\^ - y [°C]'); legend('Location','best');
xlim(XLIM_Q6);
save_v2(fig,'Q6_yerr_vs_deltaE_v2.png'); close(fig);

% =========================================================================
% Q7 - full closed loop  (recoloured lines, full time axis)
% =========================================================================
fprintf('Q7 v2...\n');

deltaE  = 1e-1;
alpha   = 1e4;
y_max7  = 55;
[A_d, B_d, C_d, L] = kalman_augmented(A, B, C, Ke, e_var, deltaE);

d_real7 = 0.50 * u_ss;
c1_p7   = c1_nom + B * d_real7;
h1_7    = h1_make(c1_p7);

r_levels = [50, 40, 60, 45];
T_step   = 600;  N_step = T_step / Ts;
N7       = N_step * numel(r_levels);
t7       = (0:N7-1)*Ts;
r7       = zeros(1,N7);
for jj = 1:numel(r_levels)
    r7((jj-1)*N_step + (1:N_step)) = r_levels(jj);
end

use_d = [true, false];
res = struct();
for cfg = 1:2
    rng(0,'twister');
    x  = nan(n,N7+1); x(:,1) = x0_amb;
    y  = nan(1,N7);   u = nan(1,N7);  dh = nan(1,N7);
    xd_hat = [zeros(n,1); 0];  Du_prev = 0;
    for k = 1:N7
        y(k)  = T1C(x(:,k));
        Dy_k  = y(k) - y_ss;
        xd_minus = A_d * xd_hat + B_d * Du_prev;
        xd_hat   = xd_minus + L * (Dy_k - C_d * xd_minus);
        Dx_hat   = xd_hat(1:n);
        d_hat    = xd_hat(n+1);
        dh(k)    = d_hat;
        Dr7 = r7(k) - y_ss;
        if use_d(cfg)
            [Dx_bar7, Du_bar7] = mpc_steady_state(A,B,C, Dr7, d_hat);
        else
            [Dx_bar7, Du_bar7] = mpc_steady_state(A,B,C, Dr7, 0);
        end
        opts = struct('formulation','dense', ...
                      'u_min',  0  - u_ss - Du_bar7, ...
                      'u_max', 100 - u_ss - Du_bar7, ...
                      'y_max', (y_max7 - y_ss) - Dr7, ...
                      'soft',  true, 'alpha', alpha);
        delta_u = mpc_solve(Dx_hat - Dx_bar7, A, B, C, H, R, opts);
        Du_k   = delta_u + Du_bar7;
        u(k)   = u_ss + Du_k;
        x(:,k+1) = h1_7(x(:,k), u(k));
        Du_prev = Du_k;
    end
    res(cfg).y = y; res(cfg).u = u; res(cfg).dh = dh;
end

fig = figure('Position',[100 100 1000 750],'Color','w');

subplot(3,1,1); hold on; grid on; box on;
title('Q7 - closed loop: reference tracking with safety cap (perturbed plant)');
plot(t7, res(1).y, 'b', 'LineWidth',1.2, 'DisplayName','y, with d\^');
plot(t7, res(2).y, 'r', 'LineWidth',1.2, 'DisplayName','y, without d\^');
stairs(t7, r7, '--', 'Color',COL_REF, 'LineWidth',1.1, 'DisplayName','Reference');
yline(y_max7,'-','Color',COL_CAP,'Label','Safety cap', ...
      'LineWidth',1.3,'HandleVisibility','off');
ylabel('y [°C]'); legend('Location','best');

subplot(3,1,2); hold on; grid on; box on;
stairs(t7, res(1).u, 'b', 'LineWidth',1.0, 'DisplayName','u, with d\^');
stairs(t7, res(2).u, 'r', 'LineWidth',1.0, 'DisplayName','u, without d\^');
yline(0,'--','Color',COL_LIMIT,'HandleVisibility','off');
yline(100,'--','Color',COL_LIMIT,'HandleVisibility','off');
ylabel('u [%]'); legend('Location','best');

subplot(3,1,3); hold on; grid on; box on;
plot(t7, res(1).dh, 'b', 'LineWidth',1.0, 'DisplayName','d\^, with d\^');
plot(t7, res(2).dh, 'r', 'LineWidth',1.0, 'DisplayName','d\^, without d\^');
yline(d_real7,'--','Color',COL_DIST, ...
      'Label',sprintf('True disturbance (%.2f)', d_real7), ...
      'LineWidth',1.3,'HandleVisibility','off');
ylim([-20, 20]);
ylabel('d\^'); xlabel('t [s]'); legend('Location','best');

save_v2(fig,'Q7_full_closed_loop_v2.png'); close(fig);

fprintf('Done. v2 figures written to %s\n', out_dir);


% =========================================================================
function [t,y,u] = run_loop_reg(A,B,C, h1, T1C, x0, x_ss, u_ss, N, Ts, H, R, opts)
    n = size(A,1);
    rng(0,'twister');
    x = nan(n,N+1); x(:,1) = x0;
    y = nan(1,N); u = nan(1,N); t = (0:N-1)*Ts;
    for k = 1:N
        y(k) = T1C(x(:,k));
        Dx_k = x(:,k) - x_ss;
        Du_k = mpc_solve(Dx_k, A, B, C, H, R, opts);
        u(k) = u_ss + Du_k;
        x(:,k+1) = h1(x(:,k), u(k));
    end
end

function [t,y,u] = run_loop_track(A,B,C, h1, T1C, x0, x_ss, u_ss, ...
                                  Dx_bar, Du_bar, N, Ts, H, R, opts)
    n = size(A,1);
    rng(0,'twister');
    x = nan(n,N+1); x(:,1) = x0;
    y = nan(1,N); u = nan(1,N); t = (0:N-1)*Ts;
    for k = 1:N
        y(k) = T1C(x(:,k));
        Dx_k = x(:,k) - x_ss;
        delta_x = Dx_k - Dx_bar;
        delta_u = mpc_solve(delta_x, A, B, C, H, R, opts);
        Du_k = delta_u + Du_bar;
        u(k) = u_ss + Du_k;
        x(:,k+1) = h1(x(:,k), u(k));
    end
end

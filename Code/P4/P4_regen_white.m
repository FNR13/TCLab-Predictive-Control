% Regenerate the P4 plots with a white background.
% Re-runs the same simulations as P4_Q2..Q5 (same rng(0) seed -> identical
% traces) and writes the PNGs straight to Imagens/P4 with white figure +
% axes background.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
common_dir = fullfile(script_dir, '..', 'Common');
addpath(common_dir);
out_dir    = fullfile(script_dir, '..', '..', '..', 'Imagens');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

% Force a light style for all figures created during this script.
set(0,'DefaultFigureColor','w');
set(0,'DefaultAxesColor','w');
set(0,'DefaultAxesXColor','k');
set(0,'DefaultAxesYColor','k');
set(0,'DefaultAxesGridColor',[0.15 0.15 0.15]);
set(0,'DefaultTextColor','k');
set(0,'DefaultAxesFontSize',11);
set(0,'DefaultLegendColor','w');
set(0,'DefaultLegendTextColor','k');
set(0,'DefaultLegendEdgeColor','k');
set(0,'DefaultColorbarColor','k');
set(0,'DefaultFigureInvertHardcopy','off');

save_white = @(fig, name) exportgraphics(fig, fullfile(out_dir,name), ...
                                         'BackgroundColor','white', ...
                                         'Resolution',150);

% ----- shared model ------------------------------------------------------
load(fullfile(common_dir,'singleheater_model.mat'), ...
     'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n     = size(A,1);
e_std = sqrt(e_var);

x_ss = [eye(n)-A; C] \ [B*u_ss; y_ss];
c1_nom = (eye(n)-A)*x_ss - B*u_ss;
c2     = y_ss - C*x_ss;

h1_make = @(c1) @(x,u) A*x + B*u + Ke*e_std*randn + c1;
T1C     = @(x)   C*x + e_std*randn + c2;

H = 20; R = 1e-2;
T = 1500;  N = T/Ts;
Dx0_amb = (eye(n)-A)\(B*(0 - u_ss));
x0_amb  = Dx0_amb + x_ss;

% =========================================================================
% Q2 - sweep H and R, unconstrained
% =========================================================================
fprintf('Regenerating Q2 plots...\n');

H_grid = [5 10 15 20 25 30];
R_grid = [1e-3, 1e-2, 1e-1, 1, 10];
nH = numel(H_grid); nR = numel(R_grid);
trace_y = cell(nH,nR); trace_u = cell(nH,nR); trace_t = cell(nH,nR);
ISE  = nan(nH,nR);  U_TV = nan(nH,nR);

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
        ISE(ih,ir)  = sum((y-y_ss).^2)*Ts;
        U_TV(ih,ir) = sum(abs(diff(u)));
    end
end

ih_plot = find(H_grid == 20);
ir_plot = find(abs(R_grid - 1e-2) < 1e-9);

% Q2 sweep R at H=20
fig = figure('Position',[100 100 900 600],'Color','w');
subplot(2,1,1); hold on; grid on; box on;
title('Q2 sweep of R at H = 20 (unconstrained)');
for ir = 1:nR
    plot(trace_t{ih_plot,ir}, trace_y{ih_plot,ir}, ...
         'LineWidth',1.0, 'DisplayName', sprintf('R = %g', R_grid(ir)));
end
yline(y_ss,'k--','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
subplot(2,1,2); hold on; grid on; box on;
for ir = 1:nR
    stairs(trace_t{ih_plot,ir}, trace_u{ih_plot,ir}, 'LineWidth',1.0);
end
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]');
save_white(fig, 'Q2_sweepR_H20.png'); close(fig);

% Q2 sweep H at R=0.01
fig = figure('Position',[100 100 900 600],'Color','w');
subplot(2,1,1); hold on; grid on; box on;
title('Q2 sweep of H at R = 0.01 (unconstrained)');
for ih = 1:nH
    plot(trace_t{ih,ir_plot}, trace_y{ih,ir_plot}, ...
         'LineWidth',1.0, 'DisplayName', sprintf('H = %d', H_grid(ih)));
end
yline(y_ss,'k--','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
subplot(2,1,2); hold on; grid on; box on;
for ih = 1:nH
    stairs(trace_t{ih,ir_plot}, trace_u{ih,ir_plot}, 'LineWidth',1.0);
end
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]');
save_white(fig, 'Q2_sweepH_R0.01.png'); close(fig);

% Q2 summary heatmap
fig = figure('Position',[100 100 900 350],'Color','w');
subplot(1,2,1);
imagesc(log10(R_grid), H_grid, log10(ISE)); colorbar;
xlabel('log_{10} R'); ylabel('H'); title('log_{10} ISE(\Delta y)');
set(gca,'YDir','normal','Color','w');
subplot(1,2,2);
imagesc(log10(R_grid), H_grid, U_TV); colorbar;
xlabel('log_{10} R'); ylabel('H'); title('TV(u)');
set(gca,'YDir','normal','Color','w');
save_white(fig, 'Q2_summary_heatmap.png'); close(fig);

% =========================================================================
% Q3 - bounded vs unbounded
% =========================================================================
fprintf('Regenerating Q3 plot...\n');

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
yline(y_ss,'k--','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
subplot(2,1,2); hold on; grid on; box on;
stairs(t_u, u_u, 'LineWidth',1.0,'DisplayName','u, unbounded');
stairs(t_b, u_b, 'LineWidth',1.2,'DisplayName','u, bounded');
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]'); legend('Location','best');
save_white(fig, 'Q3_bounded_vs_unbounded.png'); close(fig);

% =========================================================================
% Q4 - tracking nominal vs perturbed
% =========================================================================
fprintf('Regenerating Q4 plot...\n');

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
plot(t_a, y_a, 'LineWidth',1.2,'DisplayName','y (nominal)');
plot(t_b4, y_b4, 'LineWidth',1.2,'DisplayName','y (+10% c1)');
yline(r_abs,'k--','r');
yline(y_ss, 'k:','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
subplot(2,1,2); hold on; grid on; box on;
stairs(t_a, u_a, 'LineWidth',1.0,'DisplayName','u (nominal)');
stairs(t_b4, u_b4, 'LineWidth',1.0,'DisplayName','u (+10% c1)');
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]'); legend('Location','best');
save_white(fig, 'Q4_tracking_nominal_vs_perturbed.png'); close(fig);

% =========================================================================
% Q5 - hard vs soft output cap
% =========================================================================
fprintf('Regenerating Q5 plot...\n');

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
yline(r_abs5,'k--','r');
yline(y_max, 'r-' ,'y_{max}');
yline(y_ss,  'k:','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]'); legend('Location','best');
subplot(2,1,2); hold on; grid on; box on;
stairs(t_h, u_h, 'LineWidth',1.0,'DisplayName','u, hard');
stairs(t_s, u_s, 'LineWidth',1.0,'DisplayName','u, soft');
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]'); legend('Location','best');
save_white(fig, 'Q5_hard_vs_soft.png'); close(fig);

fprintf('Done. Files in %s\n', out_dir);


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

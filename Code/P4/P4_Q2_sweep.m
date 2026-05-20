% P4 Q2 - closed-loop unconstrained MPC regulator.
%
% Sweep H and R, simulate the identified TCLab model in closed loop with
% the MPC computing Du = mpc_solve(Dx, ...). No control bounds, no output
% bounds, no Kalman filter (we still use the simulator state directly).
%
% Outputs:
%   * Console table:  H, R, ISE on Dy, total |Du| energy, mean QP time.
%   * Two figure files in P4_figs/:
%       Q2_sweepR_H<H>.png  - Dy and Du for several R at the chosen H.
%       Q2_sweepH_R<R>.png  - Dy and Du for several H at the chosen R.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
common_dir = fullfile(script_dir, '..', 'Common');
addpath(common_dir);
fig_dir    = fullfile(script_dir, '..', '..', '..', 'Imagens');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% ----- model -------------------------------------------------------------
load(fullfile(common_dir,'singleheater_model.mat'), ...
     'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n     = size(A,1);
e_std = sqrt(e_var);

x_ss = [eye(n)-A; C] \ [B*u_ss; y_ss];
c1   = (eye(n)-A)*x_ss - B*u_ss;
c2   = y_ss - C*x_ss;
h1   = @(x,u) A*x + B*u + Ke*e_std*randn + c1;
T1C  = @(x)   C*x + e_std*randn + c2;

% ----- simulation common params -----------------------------------------
T = 1500;             % experiment duration [s]
N = T/Ts;             % number of samples

% Start away from equilibrium: at ambient (u = 0).
Dx0 = [eye(n)-A]\(B*(0-u_ss));   % steady-state Dx for u = 0
x0  = Dx0 + x_ss;

% ----- sweep settings ---------------------------------------------------
H_grid = [5 10 15 20 25 30];
R_grid = [1e-3, 1e-2, 1e-1, 1, 10];

% Run all combinations, store summary metrics.
nH = numel(H_grid); nR = numel(R_grid);
ISE   = nan(nH,nR);   % integral of squared Dy
U_TV  = nan(nH,nR);   % total variation of u (|Du[k]-Du[k-1]|)
SatHi = nan(nH,nR);   % fraction of samples where u > 100  (constraint *would* be hit)
SatLo = nan(nH,nR);   % fraction of samples where u <   0
Tqp   = nan(nH,nR);   % mean quadprog time [ms]

% For plotting later: keep traces at a fixed R while sweeping H, and vice
% versa. Picks are made once we see the metrics; for now just store them all.
trace_y  = cell(nH,nR);
trace_u  = cell(nH,nR);
trace_t  = cell(nH,nR);

fprintf('--- P4 Q2: closed-loop unconstrained MPC ---\n');
fprintf('%3s %7s | %10s %10s %8s %8s %9s\n', ...
        'H','R','ISE(Dy)','TV(u)','satHi%','satLo%','tQP[ms]');

rng(0,'twister');   % reproducible noise

for ih = 1:nH
    for ir = 1:nR
        H = H_grid(ih); R = R_grid(ir);

        rng(0,'twister');   % same noise sequence across designs
        x = nan(n,N+1); x(:,1) = x0;
        y = nan(1,N);  u = nan(1,N); Dy = nan(1,N); Du = nan(1,N);
        t = (0:N-1)*Ts;
        tq = zeros(1,N);

        opts = struct('formulation','dense');   % unconstrained dense

        for k = 1:N
            y(k)  = T1C(x(:,k));
            Dy(k) = y(k) - y_ss;
            Dx_k  = x(:,k) - x_ss;

            tic;
            Du_k = mpc_solve(Dx_k, A, B, C, H, R, opts);
            tq(k) = toc;

            Du(k)    = Du_k;
            u(k)     = u_ss + Du_k;
            x(:,k+1) = h1(x(:,k), u(k));
        end

        ISE(ih,ir)   = sum(Dy.^2) * Ts;
        U_TV(ih,ir)  = sum(abs(diff(Du)));
        SatHi(ih,ir) = mean(u > 100) * 100;
        SatLo(ih,ir) = mean(u <   0) * 100;
        Tqp(ih,ir)   = mean(tq) * 1e3;

        trace_y{ih,ir} = y;
        trace_u{ih,ir} = u;
        trace_t{ih,ir} = t;

        fprintf('%3d %7g | %10.2f %10.2f %8.1f %8.1f %9.2f\n', ...
                H, R, ISE(ih,ir), U_TV(ih,ir), SatHi(ih,ir), SatLo(ih,ir), Tqp(ih,ir));
    end
end

% ----- pick a representative H and R for the plots ----------------------
H_plot = 20;          % candidate horizon for Q2 figure
R_plot = 1e-2;        % candidate weight for Q2 figure
[~, ih_plot] = min(abs(H_grid - H_plot));
[~, ir_plot] = min(abs(log10(R_grid) - log10(R_plot)));
H_plot = H_grid(ih_plot);
R_plot = R_grid(ir_plot);

% ---- figure 1: sweep R at fixed H --------------------------------------
fig1 = figure('Visible','off','Position',[100 100 900 600]);
subplot(2,1,1); hold on; grid on;
title(sprintf('Q2 sweep of R at H = %d (unconstrained)', H_plot));
for ir = 1:nR
    plot(trace_t{ih_plot,ir}, trace_y{ih_plot,ir}, ...
         'LineWidth',1.0, 'DisplayName', sprintf('R = %g', R_grid(ir)));
end
yline(y_ss,'k--','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
for ir = 1:nR
    stairs(trace_t{ih_plot,ir}, trace_u{ih_plot,ir}, 'LineWidth',1.0);
end
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]');
saveas(fig1, fullfile(fig_dir, sprintf('Q2_sweepR_H%d.png', H_plot)));

% ---- figure 2: sweep H at fixed R --------------------------------------
fig2 = figure('Visible','off','Position',[100 100 900 600]);
subplot(2,1,1); hold on; grid on;
title(sprintf('Q2 sweep of H at R = %g (unconstrained)', R_plot));
for ih = 1:nH
    plot(trace_t{ih,ir_plot}, trace_y{ih,ir_plot}, ...
         'LineWidth',1.0, 'DisplayName', sprintf('H = %d', H_grid(ih)));
end
yline(y_ss,'k--','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
for ih = 1:nH
    stairs(trace_t{ih,ir_plot}, trace_u{ih,ir_plot}, 'LineWidth',1.0);
end
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]');
saveas(fig2, fullfile(fig_dir, sprintf('Q2_sweepH_R%g.png', R_plot)));

% ---- summary heatmap ---------------------------------------------------
fig3 = figure('Visible','off','Position',[100 100 900 350]);
subplot(1,2,1);
imagesc(log10(R_grid), H_grid, log10(ISE)); colorbar;
xlabel('log_{10} R'); ylabel('H'); title('log_{10} ISE(\Delta y)');
set(gca,'YDir','normal');
subplot(1,2,2);
imagesc(log10(R_grid), H_grid, U_TV); colorbar;
xlabel('log_{10} R'); ylabel('H'); title('TV(u)');
set(gca,'YDir','normal');
saveas(fig3, fullfile(fig_dir, 'Q2_summary_heatmap.png'));

fprintf('\nFigures written to %s\n', fig_dir);
fprintf('chosen H_plot = %d, R_plot = %g\n', H_plot, R_plot);

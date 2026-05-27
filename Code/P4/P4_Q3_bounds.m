% P4 Q3 - closed-loop MPC regulator with input bounds u in [0, 100]%.
%
% Uses the design fixed in Q2:  H = 20, R = 0.01.
% Compares the previously unconstrained run (Q2) with the now-bounded
% solution, on the same noise seed. We expect the unconstrained controller
% to demand u >> 100 in the first samples and the bounded one to clip at
% 100 % until the state recovers.

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

% ----- design ------------------------------------------------------------
H = 20;
R = 1e-2;

% Bounds expressed in the MPC's decision frame (increment about u_ss):
Du_min = 0   - u_ss;     % = -25
Du_max = 100 - u_ss;     % = +75

opts_unb = struct('formulation','dense');
opts_box = struct('formulation','dense', 'u_min', Du_min, 'u_max', Du_max);

% ----- simulation common params -----------------------------------------
T = 1500;
N = T/Ts;
Dx0 = (eye(n)-A)\(B*(0-u_ss));
x0  = Dx0 + x_ss;

% --- run unbounded -------------------------------------------------------
[t_u, y_u, u_u, Du_u, tq_u, ef_u] = closed_loop( ...
    A,B,C, h1, T1C, x0, x_ss, u_ss, y_ss, N, Ts, H, R, opts_unb);

% --- run bounded ---------------------------------------------------------
[t_b, y_b, u_b, Du_b, tq_b, ef_b] = closed_loop( ...
    A,B,C, h1, T1C, x0, x_ss, u_ss, y_ss, N, Ts, H, R, opts_box);

% ----- summary -----------------------------------------------------------
fprintf('--- P4 Q3: H = %d, R = %g, bounds u in [0,100] ---\n', H, R);
fprintf('                       unconstrained     bounded\n');
fprintf('peak u demanded   : %12.2f %%   %10.2f %%\n', max(u_u), max(u_b));
fprintf('min  u demanded   : %12.2f %%   %10.2f %%\n', min(u_u), min(u_b));
fprintf('samples at u_max  : %12d        %10d\n', sum(u_u >= 100), sum(u_b >= 100));
fprintf('samples at u_min  : %12d        %10d\n', sum(u_u <=   0), sum(u_b <=   0));
fprintf('ISE(Dy)           : %12.2f       %10.2f\n', sum((y_u-y_ss).^2)*Ts, sum((y_b-y_ss).^2)*Ts);
fprintf('mean tQP [ms]     : %12.2f       %10.2f\n', mean(tq_u)*1e3, mean(tq_b)*1e3);
fprintf('any exitflag != 1 : %12d        %10d\n', sum(ef_u~=1), sum(ef_b~=1));

% ----- figure ------------------------------------------------------------
fig = figure('Visible','off','Position',[100 100 950 650]);

subplot(2,1,1); hold on; grid on;
title(sprintf('Q3 - closed loop, H=%d, R=%g', H, R));
plot(t_u, y_u, 'LineWidth', 1.2, 'DisplayName','y, unbounded');
plot(t_b, y_b, 'LineWidth', 1.2, 'DisplayName','y, bounded [0,100]');
yline(y_ss,'k--','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
stairs(t_u, u_u, 'LineWidth', 1.0, 'DisplayName','u, unbounded');
stairs(t_b, u_b, 'LineWidth', 1.2, 'DisplayName','u, bounded');
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]');
legend('Location','best');

saveas(fig, fullfile(fig_dir, 'Q3_bounded_vs_unbounded.png'));
fprintf('figure written: %s\n', fullfile(fig_dir,'Q3_bounded_vs_unbounded.png'));


% =========================================================================
function [t, y, u, Du, tq, ef] = closed_loop(A,B,C, h1, T1C, x0, x_ss, ...
                                            u_ss, y_ss, N, Ts, H, R, opts)
    n = size(A,1);
    rng(0,'twister');                       % reproducible noise
    x = nan(n, N+1); x(:,1) = x0;
    y  = nan(1,N);   u  = nan(1,N);
    Du = nan(1,N);   tq = nan(1,N);   ef = nan(1,N);
    t  = (0:N-1)*Ts;
    for k = 1:N
        y(k)  = T1C(x(:,k));
        Dx_k  = x(:,k) - x_ss;

        tic;
        [Du_k, info] = mpc_solve(Dx_k, A, B, C, H, R, opts);
        tq(k) = toc;
        ef(k) = info.exitflag;

        Du(k)    = Du_k;
        u(k)     = u_ss + Du_k;
        x(:,k+1) = h1(x(:,k), u(k));
    end
end

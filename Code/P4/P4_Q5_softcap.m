% P4 Q5 - safety cap on the output: y <= y_max = 55 C.
%
% Reference r = 60 C, which is *above* the cap, so the constraint is
% guaranteed to bind. Two runs:
%   (a) hard constraint:  y_hat(i) <= y_max  enforced exactly.
%       Expected: quadprog ends up infeasible whenever the predicted
%       free-response output already crosses y_max, hence exitflag != 1
%       (or, equivalently, the QP solver returns an empty solution and
%       mpc_solve falls back to u = u_ss).
%   (b) soft constraint:  y_hat(i) <= y_max + eta_i, eta_i >= 0, penalty
%       alpha = 1e4. Always feasible; the slack absorbs the violation
%       smoothly.
%
% Mismatch turned OFF for Q5 (per guide).  H = 20, R = 1e-2.

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
H = 20; R = 1e-2;
r_abs  = 60;                          % above the cap
Dr     = r_abs - y_ss;                % +13.68
[Dx_bar, Du_bar] = mpc_steady_state(A, B, C, Dr, 0);

% Bounds (shifted)
Du_min = 0   - u_ss - Du_bar;
Du_max = 100 - u_ss - Du_bar;

% Output cap: y <= y_max  <=>  Dy <= Dy_max  <=>  δy <= Dy_max - Dr
y_max  = 55;
dy_max = (y_max - y_ss) - Dr;         % becomes NEGATIVE because r > y_max

opts_hard = struct('formulation','dense', ...
                   'u_min', Du_min, 'u_max', Du_max, ...
                   'y_max', dy_max, 'soft', false);

opts_soft = struct('formulation','dense', ...
                   'u_min', Du_min, 'u_max', Du_max, ...
                   'y_max', dy_max, 'soft', true, ...
                   'alpha', 1e4);

% ----- simulation common params -----------------------------------------
T = 1500;
N = T/Ts;
% Start ABOVE the cap to force infeasibility of the hard constraint:
% pick an equilibrium for u_hot so that y_init ~ 65 C > y_max = 55 C.
% At that initial state the plant cannot cool below y_max within the
% prediction horizon (u_min = 0, no active cooling), so the QP becomes
% infeasible.
u_hot = 70;
Dx0   = (eye(n)-A)\(B*(u_hot - u_ss));
x0    = Dx0 + x_ss;

% ----- run hard ----------------------------------------------------------
[t_h, y_h, u_h, ef_h] = closed_loop_q5( ...
    A,B,C, h1, T1C, x0, x_ss, u_ss, ...
    Dx_bar, Du_bar, N, Ts, H, R, opts_hard);

% ----- run soft ----------------------------------------------------------
[t_s, y_s, u_s, ef_s] = closed_loop_q5( ...
    A,B,C, h1, T1C, x0, x_ss, u_ss, ...
    Dx_bar, Du_bar, N, Ts, H, R, opts_soft);

% ----- summary -----------------------------------------------------------
fprintf('--- P4 Q5: cap y <= %.1f C, r = %.1f C (Dr = %.2f) ---\n', y_max, r_abs, Dr);
fprintf('                       hard          soft\n');
fprintf('exitflag == 1 always : %3s           %3s\n', tf(all(ef_h==1)), tf(all(ef_s==1)));
fprintf('# samples flag != 1  : %5d         %5d\n', sum(ef_h~=1), sum(ef_s~=1));
fprintf('unique exitflags     : ');
disp(unique(ef_h)'); fprintf('soft: '); disp(unique(ef_s)');
fprintf('max y observed       : %8.3f C   %8.3f C\n', max(y_h), max(y_s));
fprintf('mean y (last 100)    : %8.3f C   %8.3f C\n', mean(y_h(end-99:end)), mean(y_s(end-99:end)));
fprintf('mean u (last 100)    : %8.2f %%   %8.2f %%\n', mean(u_h(end-99:end)), mean(u_s(end-99:end)));

% ----- figure -----------------------------------------------------------
fig = figure('Visible','off','Position',[100 100 950 700]);

subplot(2,1,1); hold on; grid on;
title(sprintf('Q5 - reference r = %.1f C with cap y_{max} = %.1f C', r_abs, y_max));
plot(t_h, y_h, 'LineWidth', 1.2, 'DisplayName','y, hard cap');
plot(t_s, y_s, 'LineWidth', 1.2, 'DisplayName','y, soft cap');
yline(r_abs,'k--','r');
yline(y_max,'r-' ,'y_{max}');
yline(y_ss, 'k:','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
stairs(t_h, u_h, 'LineWidth', 1.0, 'DisplayName','u, hard');
stairs(t_s, u_s, 'LineWidth', 1.0, 'DisplayName','u, soft');
yline(0,'r--'); yline(100,'r--');
ylabel('u [%]'); xlabel('t [s]');
legend('Location','best');

saveas(fig, fullfile(fig_dir,'Q5_hard_vs_soft.png'));
fprintf('figure written: %s\n', fullfile(fig_dir,'Q5_hard_vs_soft.png'));


% =========================================================================
function [t, y, u, ef] = closed_loop_q5(A,B,C, h1, T1C, x0, x_ss, u_ss, ...
                                        Dx_bar, Du_bar, N, Ts, H, R, opts)
    n = size(A,1);
    rng(0,'twister');
    x = nan(n,N+1); x(:,1) = x0;
    y = nan(1,N);   u = nan(1,N);
    ef = nan(1,N);  t = (0:N-1)*Ts;
    for k = 1:N
        y(k)    = T1C(x(:,k));
        Dx_k    = x(:,k) - x_ss;
        delta_x = Dx_k - Dx_bar;
        [delta_u, info] = mpc_solve(delta_x, A, B, C, H, R, opts);
        ef(k)   = info.exitflag;
        Du_k    = delta_u + Du_bar;
        u(k)    = u_ss + Du_k;
        x(:,k+1) = h1(x(:,k), u(k));
    end
end

function s = tf(b)
    if b, s = 'yes'; else, s = 'NO'; end
end

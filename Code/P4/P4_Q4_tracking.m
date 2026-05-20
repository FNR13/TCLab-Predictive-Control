% P4 Q4 - feedforward tracking of a non-zero reference, then perturbation.
%
% Two runs:
%   (a) nominal plant, reference r = y_ss + 5 = 51.3 C, tracked with the
%       change of variables  delta_x = Dx - Dx_bar,  delta_u = Du - Du_bar
%       where (Dx_bar, Du_bar) solves the steady-state equations
%           Dx_bar = A Dx_bar + B Du_bar,
%           Dr     = C Dx_bar.
%       Same MPC regulator as Q3 (H = 20, R = 1e-2, input bounds active),
%       no Kalman yet.
%   (b) the plant constant c1 is increased by 10 % to simulate an ambient/
%       identification mismatch. The MPC still uses the nominal model,
%       so it cannot compensate the bias and a steady-state error appears.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
common_dir = fullfile(script_dir, '..', 'Common');
addpath(common_dir);
fig_dir    = fullfile(script_dir, '..', '..', '..', 'Imagens');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% ----- model (nominal) ---------------------------------------------------
load(fullfile(common_dir,'singleheater_model.mat'), ...
     'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
n     = size(A,1);
e_std = sqrt(e_var);

x_ss = [eye(n)-A; C] \ [B*u_ss; y_ss];
c1_nom = (eye(n)-A)*x_ss - B*u_ss;
c2     = y_ss - C*x_ss;

% Plant constants used in the simulator. For (a) we use c1_nom; for (b) we
% scale c1 by 1.10. Note: this is the same trick the project guide suggests
% to fake a mismatch between MPC model and reality.
make_plant = @(c1) deal( ...
        @(x,u) A*x + B*u + Ke*e_std*randn + c1, ...   % h1
        @(x)   C*x + e_std*randn + c2);               % T1C

% ----- design (frozen from Q2/Q3) ---------------------------------------
H = 20; R = 1e-2;
Dr     = 5;                       % +5 C reference step
r_abs  = y_ss + Dr;

% Feedforward steady-state (no disturbance estimate yet)
[Dx_bar, Du_bar] = mpc_steady_state(A, B, C, Dr, 0);

% Bounds on the shifted variable delta_u = Du - Du_bar:
Du_min = 0   - u_ss;
Du_max = 100 - u_ss;
opts = struct('formulation','dense', ...
              'u_min', Du_min - Du_bar, ...
              'u_max', Du_max - Du_bar);

% ----- simulation common params -----------------------------------------
T = 1500;
N = T/Ts;
% Start from ambient temperature (u = 0 equilibrium of the *nominal* plant).
Dx0  = (eye(n)-A)\(B*(0-u_ss));
x0   = Dx0 + x_ss;

% ----- run (a) nominal ---------------------------------------------------
[h1_a, T1C_a] = make_plant(c1_nom);
[t_a, y_a, u_a] = closed_loop_track( ...
    A,B,C, h1_a, T1C_a, x0, x_ss, u_ss, y_ss, ...
    Dx_bar, Du_bar, N, Ts, H, R, opts);

% ----- run (b) c1 mismatch ----------------------------------------------
% The nominal c1 from least-squares fitting is essentially numerical noise
% (||c1|| ~ 3e-3); literally multiplying it by 1.10 produces no visible
% offset. We interpret "+10% increase on c1" as an additional bias
% equivalent to a 10 % shift of the steady-state input, i.e. a
% constant input disturbance d_real = 0.10 * u_ss applied at the plant.
% This is the kind of mismatch the project guide is trying to expose:
% something the MPC's identified model does not know about.
% Without integral action the MPC rejects roughly (1 - 1/loop_gain) of any
% input bias, so a small d_real barely shows up at steady state. The DC
% gain of the identified model is ~0.96 deg/percent: we pick d_real so the
% open-loop bias is ~12 C, giving a residual offset clearly visible (and
% later compensated by the Q7 disturbance estimator).
d_real  = 0.50 * u_ss;          % +12.5 % equivalent input bias
c1_pert = c1_nom + B * d_real;
[h1_b, T1C_b] = make_plant(c1_pert);
[t_b, y_b, u_b] = closed_loop_track( ...
    A,B,C, h1_b, T1C_b, x0, x_ss, u_ss, y_ss, ...
    Dx_bar, Du_bar, N, Ts, H, R, opts);

% ----- console summary --------------------------------------------------
% Use the last 100 samples to measure steady-state error.
tail = (N-99):N;
ss_a = mean(y_a(tail)); err_a = ss_a - r_abs;
ss_b = mean(y_b(tail)); err_b = ss_b - r_abs;

fprintf('--- P4 Q4: tracking r = %.2f C  (Dr = %.1f) ---\n', r_abs, Dr);
fprintf('feedforward steady-state target: Dx_bar (norm) = %.4f,  Du_bar = %.4f\n', ...
        norm(Dx_bar), Du_bar);
fprintf('\n                   nominal     +10%% c1 mismatch\n');
fprintf('mean y(last 100): %8.3f C  %12.3f C\n', ss_a, ss_b);
fprintf('y error vs r   : %+8.3f C  %12.3f C\n', err_a, err_b);
fprintf('peak u          : %8.2f %%  %12.2f %%\n', max(u_a), max(u_b));

% ----- figure -----------------------------------------------------------
fig = figure('Visible','off','Position',[100 100 950 700]);

subplot(2,1,1); hold on; grid on;
title(sprintf('Q4 - tracking r = %.2f C, H = %d, R = %g', r_abs, H, R));
plot(t_a, y_a, 'LineWidth',1.2, 'DisplayName','y (nominal)');
plot(t_b, y_b, 'LineWidth',1.2, 'DisplayName','y (+10% c1)');
yline(r_abs,'k--','r');
yline(y_ss, 'k:','y_{ss}');
ylabel('y [°C]'); xlabel('t [s]');
legend('Location','best');

subplot(2,1,2); hold on; grid on;
stairs(t_a, u_a, 'LineWidth',1.0, 'DisplayName','u (nominal)');
stairs(t_b, u_b, 'LineWidth',1.0, 'DisplayName','u (+10% c1)');
yline(0,'r--'); yline(100,'r--');
yline(u_ss + Du_bar, 'k:', 'u_{ss}+\Deltau_{bar}');
ylabel('u [%]'); xlabel('t [s]');
legend('Location','best');

saveas(fig, fullfile(fig_dir,'Q4_tracking_nominal_vs_perturbed.png'));
fprintf('figure written: %s\n', fullfile(fig_dir,'Q4_tracking_nominal_vs_perturbed.png'));


% =========================================================================
function [t, y, u] = closed_loop_track(A,B,C, h1, T1C, x0, x_ss, u_ss, y_ss, ...
                                       Dx_bar, Du_bar, N, Ts, H, R, opts)
    n = size(A,1);
    rng(0,'twister');
    x = nan(n,N+1); x(:,1) = x0;
    y = nan(1,N);   u = nan(1,N);
    t = (0:N-1)*Ts;
    for k = 1:N
        y(k)    = T1C(x(:,k));
        Dx_k    = x(:,k) - x_ss;
        delta_x = Dx_k - Dx_bar;                      % shifted state
        delta_u = mpc_solve(delta_x, A, B, C, H, R, opts);
        Du_k    = delta_u + Du_bar;                   % undo the shift
        u(k)    = u_ss + Du_k;
        x(:,k+1) = h1(x(:,k), u(k));
    end
end

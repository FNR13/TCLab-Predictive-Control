% P4 Q1 - sanity check: mpc_solve must reproduce the analytical P2 gain.
%
% For the scalar plants of P2 (A in {1.2, 0.8}, B=C=1, Q=1), the optimal
% receding-horizon control is the linear state feedback
%       u0_analytic = -K_RH * x0,        K_RH = e1 * (W'W + R I)^-1 * W' * Pi.
%
% The script
%   * computes K_RH analytically for a grid of (A, R, H),
%   * calls mpc_solve in dense and sparse mode for x0 = 1,
%   * compares u0 to -K_RH and reports the worst-case absolute error.
%
% No figures: the test is numerical and the verdict is printed.

clear; close all; clc;

% Locate the shared Common/ folder (mpc_solve and helpers live there).
script_dir = fileparts(mfilename('fullpath'));
common_dir = fullfile(script_dir, '..', 'Common');
addpath(common_dir);

A_list = [1.2, 0.8];
R_list = [0.1, 1, 10, 100];
H_list = [1, 5, 10, 20, 30];
B = 1; C = 1;
x0 = 1;

fprintf('--- P4 Q1: mpc_solve vs P2 analytic K_RH ---\n');
fprintf('%4s %5s %4s | %-13s %-13s %-13s %-13s\n', ...
        'A','R','H','K_RH (ana.)','u0 dense','u0 sparse','max |err|');

worst = 0;
for A = A_list
    for R = R_list
        for H = H_list
            % Analytic K_RH (same construction as in P2)
            W  = zeros(H); Pi = zeros(H,1);
            for i = 1:H
                Pi(i) = C * A^i;
                for j = 1:i
                    W(i,j) = C * A^(i-j) * B;
                end
            end
            M     = W.'*W + R*eye(H);
            K_RH  = [1, zeros(1,H-1)] * (M \ (W.' * Pi));   % scalar
            u_ana = -K_RH * x0;

            % mpc_solve - dense
            opts_d = struct('formulation','dense');
            u_d = mpc_solve(x0, A, B, C, H, R, opts_d);

            % mpc_solve - sparse
            opts_s = struct('formulation','sparse');
            u_s = mpc_solve(x0, A, B, C, H, R, opts_s);

            err = max(abs([u_d - u_ana, u_s - u_ana]));
            worst = max(worst, err);

            fprintf('%4.1f %5g %4d | %+13.6e %+13.6e %+13.6e %13.3e\n', ...
                    A, R, H, K_RH, u_d, u_s, err);
        end
    end
end

fprintf('\nWorst-case |u0 - u_analytic| across all combinations: %.3e\n', worst);
if worst < 1e-6
    fprintf('PASS: mpc_solve matches the analytical P2 result.\n');
else
    fprintf('FAIL: discrepancy above 1e-6 tolerance.\n');
end

function [Dx_bar, Du_bar] = mpc_steady_state(A, B, C, Dr, d_hat)
% MPC_STEADY_STATE  Feedforward target for tracking with input disturbance.
%
% Solves
%     Dx_bar = A Dx_bar + B (Du_bar + d_hat)
%     Dr     = C Dx_bar
% for (Dx_bar, Du_bar). Used in P4 Q4 (with d_hat = 0) and Q7 (with the
% Kalman disturbance estimate).
%
% Inputs
%   A,B,C  state-space matrices of the incremental model.
%   Dr     scalar reference increment (Dr = r - y_ss).
%   d_hat  scalar input-disturbance estimate (use 0 if not estimated).
%
% Outputs
%   Dx_bar  n x 1 steady-state increment of the state.
%   Du_bar  scalar steady-state increment of the control.

if nargin < 5 || isempty(d_hat), d_hat = 0; end

n = size(A,1);
M  = [eye(n) - A,  -B;
            C   ,   0];
rhs = [B * d_hat;  Dr];

sol    = M \ rhs;
Dx_bar = sol(1:n);
Du_bar = sol(n+1);
end

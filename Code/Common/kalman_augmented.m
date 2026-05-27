function [A_d, B_d, C_d, L] = kalman_augmented(A, B, C, Ke, e_var, delta_E)
% KALMAN_AUGMENTED  Build the disturbance-augmented Kalman filter for P4 Q6.
%
% The augmented state is
%     x_d = [Delta_x ; d]
% with dynamics
%     x_d(k+1) = A_d x_d(k) + B_d Delta_u(k) + noise
%     Delta_y  = C_d x_d(k) + e(k)
% and the augmented process-noise covariance
%     Q_E_d = blkdiag( K_e * e_var * K_e' ,  delta_E ).
%
% Inputs
%   A,B,C    incremental-model matrices (n x n, n x 1, 1 x n).
%   Ke       innovation gain from ssest (n x 1).
%   e_var    measurement-noise variance from ssest (scalar).
%   delta_E  tuning variance for the disturbance state (scalar > 0).
%
% Outputs
%   A_d, B_d, C_d  augmented state-space matrices (size n+1).
%   L              steady-state Kalman gain ((n+1) x 1).
%
% Use the returned (A_d, B_d, C_d, L) in the simulator/controller loop:
%
%   % prediction
%   x_d_minus = A_d * x_d_prev + B_d * Du_prev;
%   % correction with the latest measurement Dy(k)
%   x_d       = x_d_minus + L * ( Dy(k) - C_d * x_d_minus );
%
%   Dx_hat = x_d(1:n);
%   d_hat  = x_d(n+1);

n   = size(A,1);
A_d = [A,  B;
       zeros(1,n),  1];
B_d = [B; 0];
C_d = [C,  0];

Q_E   = (Ke * Ke.') * e_var;             % n x n
Q_E_d = blkdiag(Q_E, delta_E);           % (n+1) x (n+1)
R_E   = e_var;                           % scalar

L = dlqe(A_d, eye(n+1), C_d, Q_E_d, R_E);
end

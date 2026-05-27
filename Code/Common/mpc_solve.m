function [u0, info] = mpc_solve(x0, A, B, C, H, R, opts)
% MPC_SOLVE  One step of a finite-horizon MPC for a SISO linear system.
%
% Solves
%     min sum_{i=0..H-1} y_hat(i+1)^2 + R * u_hat(i)^2
%     s.t.  x_hat(i+1) = A x_hat(i) + B u_hat(i)
%           y_hat(i)   = C x_hat(i)
%           u_min <= u_hat(i) <= u_max         (optional, P4 Q3)
%           y_hat(i) <= y_max + eta_i          (optional, P4 Q5 soft)
%           eta_i >= 0
% with x_hat(0) = x0, and returns only the first control u0 = u_hat(0).
%
% The function supports two equivalent formulations selected through
% opts.formulation: 'dense' (decision vector = U) or 'sparse'
% (decision vector = [X; U]). Soft output constraints append eta to the
% decision vector in either case.
%
% Inputs
%   x0    n x 1 initial (incremental) state, x0 = Delta_x(k)  [or delta_x(k)
%         after the change of variables for tracking, see Section 5 of
%         P4_derivations.md]
%   A,B,C state-space matrices of the incremental model (Delta dynamics).
%         B is n x 1, C is 1 x n.
%   H     prediction horizon (number of control moves)
%   R     scalar control weight (R > 0)
%   opts  struct with optional fields:
%           .formulation  'dense' (default) or 'sparse'
%           .u_min        scalar, lower bound on the optimisation input
%                         (already in the same frame as the decision var;
%                         e.g. for Q3 pass Delta_u_min = u_min - u_ss;
%                         for tracking pass Delta_u_min - Delta_u_bar)
%           .u_max        scalar, upper bound
%           .y_max        scalar, upper bound on the optimisation output
%                         (same frame as the decision: typically
%                         Delta_y_max - Delta_r in the tracking case)
%           .soft         logical, true to relax y_max into a soft
%                         constraint with slack eta_i (default false)
%           .alpha        scalar, soft-constraint penalty (default 1e4)
%           .verbose      logical, print quadprog exitflag info
%
% Outputs
%   u0    scalar, first control move to apply to the plant
%         (in the same frame as the optimisation variable)
%   info  struct with fields:
%           .U            full optimal control sequence (H x 1)
%           .exitflag     quadprog exit flag
%           .formulation  which path was used
%           .solve_time   seconds spent in quadprog
%
% Notes
%   * For Q2 (unconstrained regulator) leave opts empty (or omit it).
%   * For Q3 (input bounds) set opts.u_min / opts.u_max.
%   * For Q4 (feedforward tracking) the caller is responsible for the
%     change of variables: pass x0 = Delta_x - Delta_x_bar and shift the
%     bounds by -Delta_u_bar / -Delta_r before calling.
%   * For Q5 (soft output bound) set opts.y_max and opts.soft = true.
%
% See also P4_derivations.md (sections 3, 4, 5).

% ---------- argument handling ----------------------------------------------
if nargin < 7 || isempty(opts), opts = struct(); end
if ~isfield(opts,'formulation'), opts.formulation = 'dense'; end
if ~isfield(opts,'soft'),        opts.soft        = false;   end
if ~isfield(opts,'alpha'),       opts.alpha       = 1e4;     end
if ~isfield(opts,'verbose'),     opts.verbose     = false;   end

has_ubox  = isfield(opts,'u_min') || isfield(opts,'u_max');
has_ybnd  = isfield(opts,'y_max');

n = size(A,1);
assert(size(B,1) == n && size(B,2) == 1, 'B must be n x 1');
assert(size(C,1) == 1 && size(C,2) == n, 'C must be 1 x n');
assert(numel(x0) == n, 'x0 must be n x 1');
x0 = x0(:);

% ---------- dispatch -------------------------------------------------------
t0 = tic;
switch lower(opts.formulation)
    case 'dense'
        [U, exitflag] = solve_dense(x0, A, B, C, H, R, opts, has_ubox, has_ybnd);
    case 'sparse'
        [U, exitflag] = solve_sparse(x0, A, B, C, H, R, opts, has_ubox, has_ybnd);
    otherwise
        error('mpc_solve:badFormulation', ...
              'opts.formulation must be ''dense'' or ''sparse''.');
end
solve_time = toc(t0);

if opts.verbose && exitflag ~= 1
    fprintf('[mpc_solve] quadprog exitflag = %d\n', exitflag);
end

u0   = U(1);
info = struct('U', U, ...
              'exitflag', exitflag, ...
              'formulation', opts.formulation, ...
              'solve_time', solve_time);
end


% ========================================================================
% DENSE FORMULATION  (decision vector: z = U [;eta])
% ========================================================================
function [U, exitflag] = solve_dense(x0, A, B, C, H, R, opts, has_ubox, has_ybnd)

[W, Pi] = build_W_Pi(A, B, C, H);            % H x H and H x n

% ---- cost on U only:  J = U' (W'W + R I) U + 2 x0' Pi' W U  (+ const) ----
F_U = 2 * (W.'*W + R*eye(H));
f_U = 2 * (W.' * Pi * x0);

% ---- assemble (with optional eta block) -----------------------------------
soft = opts.soft && has_ybnd;

if soft
    F = blkdiag(F_U, 2*opts.alpha*eye(H));   % factor 2 to match quadprog 1/2 z'Fz
    f = [f_U; zeros(H,1)];
else
    F = F_U;
    f = f_U;
end
F = (F + F.')/2;                             % enforce numerical symmetry

% ---- inequality constraints (hard or soft output cap) ---------------------
A_ineq = [];
b_ineq = [];
if has_ybnd
    % y_hat <= y_max  =>  W U <= y_max - Pi x0  (+ eta if soft)
    g_tilde = opts.y_max * ones(H,1) - Pi * x0;
    if soft
        A_ineq = [W,  -eye(H)];
        b_ineq = g_tilde;
    else
        A_ineq = W;
        b_ineq = g_tilde;
    end
end

% ---- bounds ---------------------------------------------------------------
if has_ubox
    u_min = getfield_default(opts,'u_min',-inf);
    u_max = getfield_default(opts,'u_max', inf);
    lb_U = u_min * ones(H,1);
    ub_U = u_max * ones(H,1);
else
    lb_U = -inf(H,1);
    ub_U =  inf(H,1);
end
if soft
    lb = [lb_U;  zeros(H,1)];                % eta >= 0
    ub = [ub_U;   inf(H,1)];
else
    lb = lb_U;
    ub = ub_U;
end

% ---- solve ----------------------------------------------------------------
qp_opts = optimoptions('quadprog','Display','off');
[z, ~, exitflag] = quadprog(F, f, A_ineq, b_ineq, [], [], lb, ub, [], qp_opts);

if isempty(z)
    % infeasible / failed; fall back to a zero control move so the loop
    % keeps running. The caller can inspect info.exitflag.
    U = zeros(H,1);
else
    U = z(1:H);
end
end


% ========================================================================
% SPARSE FORMULATION  (decision vector: z = [X; U] [;eta])
% ========================================================================
function [U, exitflag] = solve_sparse(x0, A, B, C, H, R, opts, has_ubox, has_ybnd)

n = size(A,1);
nX = (H+1)*n;       % rows of X block in z
nU = H;             % rows of U block in z

% ---- dynamics as equality constraint -------------------------------------
%   X = A_tilde X + B_tilde U + E x0   <=>   (I - A_tilde) X - B_tilde U = E x0
A_tilde = build_A_tilde(A, H, n);
B_tilde = build_B_tilde(B, H, n);
E       = [eye(n); zeros(H*n, n)];

A_eq = [eye(nX) - A_tilde,  -B_tilde];
b_eq = E * x0;

% ---- cost ----------------------------------------------------------------
%   J = X' Q_tilde X + U' R_tilde U     with Q_tilde = blkdiag(0, Q, ..., Q)
Q = C.' * C;
Q_tilde = blkdiag(zeros(n), kron(eye(H), Q));   % (H+1)n x (H+1)n
R_tilde = R * eye(H);

% ---- optional eta block --------------------------------------------------
soft = opts.soft && has_ybnd;
if soft
    F = 2 * blkdiag(Q_tilde, R_tilde, opts.alpha*eye(H));
    f = zeros(nX + nU + H, 1);
    A_eq = [A_eq, zeros(nX, H)];
else
    F = 2 * blkdiag(Q_tilde, R_tilde);
    f = zeros(nX + nU, 1);
end
F = (F + F.')/2;

% ---- inequality constraints (output) -------------------------------------
A_ineq = [];
b_ineq = [];
if has_ybnd
    C_tilde = build_C_tilde(C, H, n);     % H x (H+1)n
    g_tilde = opts.y_max * ones(H,1);
    if soft
        A_ineq = [C_tilde, zeros(H, nU), -eye(H)];
    else
        A_ineq = [C_tilde, zeros(H, nU)];
    end
    b_ineq = g_tilde;
end

% ---- bounds --------------------------------------------------------------
if has_ubox
    u_min = getfield_default(opts,'u_min',-inf);
    u_max = getfield_default(opts,'u_max', inf);
    lb_U = u_min * ones(H,1);
    ub_U = u_max * ones(H,1);
else
    lb_U = -inf(H,1);
    ub_U =  inf(H,1);
end
if soft
    lb = [-inf(nX,1); lb_U; zeros(H,1)];
    ub = [ inf(nX,1); ub_U;  inf(H,1)];
else
    lb = [-inf(nX,1); lb_U];
    ub = [ inf(nX,1); ub_U];
end

% ---- solve ---------------------------------------------------------------
qp_opts = optimoptions('quadprog','Display','off');
[z, ~, exitflag] = quadprog(F, f, A_ineq, b_ineq, A_eq, b_eq, lb, ub, [], qp_opts);

if isempty(z)
    U = zeros(H,1);
else
    U = z(nX + (1:nU));
end
end


% ========================================================================
% Building blocks
% ========================================================================
function [W, Pi] = build_W_Pi(A, B, C, H)
% W is the Toeplitz prediction matrix mapping U -> Y.
% Pi maps the initial state x0 -> Y (free response).
n  = size(A,1);
W  = zeros(H,H);
Pi = zeros(H,n);

A_pow = eye(n);                  % A^0
% column j contains the impulse response shifted by j-1 samples.
% Pre-compute CA^i B for i = 0..H-1, and CA^(i+1) for Pi.
CAB = zeros(1,H);
for i = 0:H-1
    CAB(i+1) = C * A_pow * B;     % C A^i B
    A_pow    = A * A_pow;         % becomes A^(i+1) after the line
    Pi(i+1,:) = C * A_pow;        % C A^(i+1)
end
% Fill W as a lower-triangular Toeplitz: W(i,j) = C A^(i-j) B for j <= i.
for i = 1:H
    for j = 1:i
        W(i,j) = CAB(i - j + 1);
    end
end
end


function A_tilde = build_A_tilde(A, H, n)
% Block-lower-shift with A on the sub-diagonal; (H+1) block rows/cols.
A_tilde = zeros((H+1)*n);
for k = 1:H
    rows = k*n + (1:n);          % block row k+1 (i.e. x_hat(k))
    cols = (k-1)*n + (1:n);      % block col k   (i.e. x_hat(k-1))
    A_tilde(rows, cols) = A;
end
end


function B_tilde = build_B_tilde(B, H, n)
% (H+1) block rows, H block cols of width 1; B on the (k+1,k) block.
B_tilde = zeros((H+1)*n, H);
for k = 1:H
    rows = k*n + (1:n);          % x_hat(k) row
    B_tilde(rows, k) = B;
end
end


function C_tilde = build_C_tilde(C, H, n)
% Selects rows x_hat(1)..x_hat(H) of X and applies C.
C_tilde = zeros(H, (H+1)*n);
for k = 1:H
    cols = k*n + (1:n);
    C_tilde(k, cols) = C;
end
end


function v = getfield_default(s, fname, default)
if isfield(s, fname) && ~isempty(s.(fname))
    v = s.(fname);
else
    v = default;
end
end

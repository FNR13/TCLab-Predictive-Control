% P1 - Constrained vs unconstrained optimisation of the Rosenbrock function.
% Reproduces Figures 1 and 2 of the project report.
%
% f(x1, x2) = 100*(x2 - x1^2)^2 + (1 - x1)^2

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, '..', '..', '..', 'Imagens');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

set(0,'DefaultFigureColor','w','DefaultAxesColor','w', ...
      'DefaultAxesXColor','k','DefaultAxesYColor','k', ...
      'DefaultAxesGridColor',[0.15 0.15 0.15],'DefaultTextColor','k', ...
      'DefaultAxesFontSize',12,'DefaultLegendColor','w', ...
      'DefaultLegendTextColor','k','DefaultLegendEdgeColor','k', ...
      'DefaultFigureInvertHardcopy','off');
save_white = @(fig, name) exportgraphics(fig, fullfile(out_dir,name), ...
                                         'BackgroundColor','white','Resolution',150);

% Initial estimate
x0 = [-1; 1];

% Unconstrained min
opts_unc = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
xopt = fminunc(@RosenbrockFunction, x0, opts_unc);
fopt = RosenbrockFunction(xopt);
fprintf('Unconstrained minimum:  x* = [%.4f, %.4f]   f(x*) = %.4e\n', ...
        xopt(1), xopt(2), fopt);

% Constrained min:  x1 <= 0.5
A_ineq = [1 0];
b_ineq = 0.5;
opts_con = optimoptions('fmincon','Algorithm','sqp','Display','off');
xoptc = fmincon(@RosenbrockFunction, x0, A_ineq, b_ineq, [], [], [], [], [], opts_con);
foptc = RosenbrockFunction(xoptc);
fprintf('Constrained minimum:    x_c* = [%.4f, %.4f]   f(x_c*) = %.4e\n', ...
        xoptc(1), xoptc(2), foptc);

% Grid for plots: same window as in the report (x1 in [-1.5, 1.5], x2 in [-0.5, 2])
N = 200;
x1v = linspace(-1.5, 1.5, N);
x2v = linspace(-0.5, 2.0, N);
[X1, X2] = meshgrid(x1v, x2v);
F = 100*(X2 - X1.^2).^2 + (1 - X1).^2;

% ======================================================================
% Figure 1 - 3D surface with the constraint plane (matches report fig. 1)
% ======================================================================
fig = figure('Position',[100 100 800 600],'Color','w');
surf(X1, X2, F, 'EdgeColor','none');
shading interp; hold on;
% Constraint plane x1 = 0.5 (red, semi-transparent)
x1p = [0.5 0.5 0.5 0.5];
x2p = [-0.5 2.0 2.0 -0.5];
fp  = [0 0 500 500];
patch(x1p, x2p, fp, [1 0.4 0.4], 'FaceAlpha', 0.35, 'EdgeColor','none');
% Match the report's viewing angle (camera looks from x1 negative, x2 positive)
view(-30, 30);
zlim([0, 500]);
xlim([-1.5, 1.5]); ylim([-0.5, 2.0]);
xlabel('x_1'); ylabel('x_2'); zlabel('f(x)');
title('Rosenbrock function with constraint plane x_1 = 0.5');
% Mark initial estimate, unconstrained and constrained minima with red markers
plot3(x0(1),    x0(2),    RosenbrockFunction(x0),    'ro', 'MarkerSize',8, 'LineWidth',1.5);
plot3(xopt(1),  xopt(2),  RosenbrockFunction(xopt),  'rx', 'MarkerSize',10,'LineWidth',2);
plot3(xoptc(1), xoptc(2), RosenbrockFunction(xoptc), 'r*', 'MarkerSize',10,'LineWidth',2);
save_white(fig, 'P1_Rosenbrock_surf.png'); close(fig);

% ======================================================================
% Figure 2 - contour plot (matches report fig. 2)
% ======================================================================
% Report uses ~10 level curves and a colorbar that goes up to ~600.
fig = figure('Position',[100 100 800 600],'Color','w');
contour(X1, X2, F, 10, 'LineWidth', 1.2); colorbar; hold on; grid on; box on;
% Constraint boundary x1 = 0.5 (solid black vertical line)
xline(0.5, 'k', 'LineWidth', 1.5);
% Markers
plot(x0(1),    x0(2),    'ro', 'MarkerSize', 10, 'LineWidth', 1.5);  % initial
plot(xopt(1),  xopt(2),  'rx', 'MarkerSize', 12, 'LineWidth', 2);    % unconstrained
plot(xoptc(1), xoptc(2), 'r*', 'MarkerSize', 12, 'LineWidth', 2);    % constrained
xlim([-1.5, 1.5]); ylim([-0.5, 2.0]);
xlabel('x_1'); ylabel('x_2');
save_white(fig, 'P1_Rosenbrock_contour.png'); close(fig);

fprintf('Figures written to %s\n', out_dir);

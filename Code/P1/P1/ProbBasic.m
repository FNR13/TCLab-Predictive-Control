% Illustrates the use of the Matlab solver fminunc to find the minimum of
% a quadratic function.
%
% Matlab sw required: optimization toolbox
% Functions called
%    fminunc - Matlab function (optimization toolbox) for uncionstrained
%        minimization
%    BasicFunction.m - user defined function that defines the function to be
%        minimized
%
% IST, MEEC, Distributed Predictive Control and Estimation
% J. Miranda Lemos, 2022
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% =========================================================
%% P1 - modifided for laboratory work 05/2026
%%
%% Group 17
%% Miguel Brás Simões - 99162
%% Henrique Miranda Póvoa - -102109
%% David de Jesus Alpalhão Marques - 100145
%% Ricardo Nascimento Francisco - 103093
%% =========================================================



% Plots function level curves

% Range of independent variables to consider 
x1min=-1.5;
x1max=1.5;
x2min=-0.5;
x2max=2;

% Number of intervals in the mesh grid
N1=100;
N2=100;

xv1 = linspace(x1min,x1max,N1);
xv2 = linspace(x2min,x2max,N2);
[xx1,xx2] = meshgrid(xv1,xv2);

% Computes the function at the different points of the mesh grid
for ii=1:N1
    for jj=1:N2
        x=[xx1(ii,jj); xx2(ii,jj)];
        ff(ii,jj)=BasicFunction(x);
    end
end

% Plots the level curves using the Matlab function contoutr
Nlevel=10;  % Number of level curves in the contour plot
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure(1), contour(xv1,xv2,ff,Nlevel,LW,1.2), colorbar
axis([x1min x1max x2min x2max]), axis square
hold on

%--------------------------------------------------------------------------
% Compute the minimum

% Initial estimate of the minimum
%x0=[0.2; 0.5];
x0=[-1; 1];
% Define the options to be used with the fminunc solver:
% The quasi-newton algorithm is used because it does not require the
% gradient; the default algorithm requires the gradient
options = optimoptions('fminunc','Algorithm','quasi-newton');

% Uses the solver fminunc to compute the minimum of the function defined in
% the Matlab function defined in the file BasicFunction.m
xopt=fminunc(@BasicFunction,x0,options)

%--------------------------------------------------------------------------
% Computes the constrained minimum associated to the constraint 
% x(1)+x(2) <= 1
%---Restrições problema 0------
%Ac=[1 1];
%Bc=1;
%------------------------------
%---Restrições problema 1------
% x(1) <= 0.5
Ac=[1 0];
Bc=0.5;
%------------------------------

xoptconstr=fmincon(@BasicFunction,x0,Ac,Bc)

%--------------------------------------------------------------------------
% Plots the initial point as a red circle
gg=plot(x0(1),x0(2),'or');
set(gg,'Linewidth',1.5);

% Plots the final estimate of the unconstrained minimum as a red cross
gg=plot(xopt(1),xopt(2),'xr');
set(gg,'Linewidth',1.5);

% Plots the final estimate of the constrained minimum as a red star
gg=plot(xoptconstr(1),xoptconstr(2),'*r');
set(gg,'Linewidth',1.5);

%plots the constraint boundary

%z1c=[0:0.1:1];
%z2c=[1:-0.1:0];
z2c = linspace(x2min,x2max,100);
z1c = 0.5*ones(size(z2c));

gg=plot(z1c,z2c,'k');
set(gg,'Linewidth',1.5);

% Identifies axis
gg=xlabel('x_1');
set(gg,'FontSize',14);

gg=ylabel('x_2');
set(gg,'FontSize',14);

hold off

%--------------------------------------------------------------------------
% Plots the 3D view of the function

figure(2)
surf(xx1,xx2,ff,'EdgeColor','k');
hold on
shading interp

% --- Single vertical constraint plane: x1 = 0.5 ---
zmin = min(ff(:));
zmax = max(ff(:));

[X2p,Zp] = meshgrid([x2min x2max],[zmin zmax]);
X1p = 0.5*ones(size(X2p));

surf(X1p,X2p,Zp, ...
    'FaceColor','red', ...
    'FaceAlpha',0.35, ...
    'EdgeColor','none');

% --- Compute function values at important points ---
f0 = BasicFunction(x0);
fopt = BasicFunction(xopt);
foptconstr = BasicFunction(xoptconstr);

% --- Mark the points on the 3D plot ---
plot3(xopt(1),xopt(2),fopt,'xr','LineWidth',2,'MarkerSize',10);          % unconstrained optimum
plot3(xoptconstr(1),xoptconstr(2),foptconstr,'*r','LineWidth',2,'MarkerSize',10); % constrained optimum
plot3(x0(1),x0(2),f0,'or','LineWidth',1.5,'MarkerSize',8);               % initial point

% Optional: dashed vertical lines to help visualize the points
plot3([xopt(1) xopt(1)],[xopt(2) xopt(2)],[zmin fopt],'--r','LineWidth',1);
plot3([xoptconstr(1) xoptconstr(1)],[xoptconstr(2) xoptconstr(2)],[zmin foptconstr],'--r','LineWidth',1);

% Identifies axis
gg = xlabel('x_1');
set(gg,'FontSize',14);

gg = ylabel('x_2');
set(gg,'FontSize',14);

gg = zlabel('f(x)');
set(gg,'FontSize',14);

title('Rosenbrock function with constraint plane x_1 = 0.5')
view(45,30)
grid on
hold off
%--------------------------------------------------------------------------
% End of File

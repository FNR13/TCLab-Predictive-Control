function f = RosenbrockFunction(x)
% Rosenbrock function in 2D:  f(x) = 100*(x2 - x1^2)^2 + (1 - x1)^2.
%
% Used in P1 to illustrate unconstrained (fminunc) and constrained
% (fmincon) optimisation. The unconstrained minimum is (1,1) with f = 0;
% with the constraint x1 <= 0.5 the optimum is (0.5, 0.25) with f = 0.25.

f = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
end

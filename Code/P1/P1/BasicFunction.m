function f = BasicFunction(x)
%   This function computes a quadratic function in 2 dimensions at a point
%   The argument x of the function is a 2-dimensional vector;
%   The function computes the value of a quadratic function at x.
%
%   IST - MEEC - Distributed Predictive Control and Estimation
%   J. Miranda Lemos, 2022
%--------------------------------------------------------------------------

%f = (x(1)-1)^2 + (x(2)-1)^2;   %função ponto 0
f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;    %função ponto 1
end
%--------------------------------------------------------------------------
% End of file



% Loadf.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the local 1D load vector for FEM.

function f = Loadf(x)

% ===Input = x, mesh nodes at which to evaluate the load function f
% ===Output = f, value of load f at x

f = (pi^2)*sin(pi*x);

end
% Exact.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the exact solution for the simple FEM method.

function y = Exact(x)

% ===Input = x, mesh nodes at which to evaluate the exact solution y
% ===Output = y, the exact solution at x

y = sin(pi*x);

end
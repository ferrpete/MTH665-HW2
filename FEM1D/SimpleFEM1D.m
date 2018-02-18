% SimpleFEM1D.m
% Peter Ferrero, Oregon State University, MTH655, 1/31/2018
% A simple FEM 1D method to solve Problem 4 of Homework 1 for MTH 655.

function [FemSol, x] = SimpleFEM1D(N)

a = 0; % left endpoint
b = 1; % right endpoint
h = (b-a)/N; % uniform mesh size
x = a:h:b; % mesh nodes

A = GStiff(x); % Global Stiffness matrix
F = GLoad(x); % Load vector

FemSol = zeros(N+1, 1); % Initialize the FEM solution
FemSol(2:N) = A(2:N,2:N)\F(2:N); % Solve the linear system
                                 % for interior nodes

end
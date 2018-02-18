% SimpleFEM1DRef.m
% Peter Ferrero, Oregon State University, MTH655, 1/31/2018
% A simple 1D reference element method to solve Problem 4 of Homework 2 for MTH 655.

function [FemSol, x] = SimpleFEM1DRef(N)

a = 0; % left endpoint
b = 1; % right endpoint
h = (b-a)/(2*N); % uniform mesh size
x = a:h:b; % mesh nodes

A = GStiffRef(x); % Global Stiffness matrix
F = GLoadRef(x); % Load vector

FemSol = zeros(2*N+1, 1); % Initialize the FEM solution
FemSol(2:end-1) = A(2:end-1,2:end-1)\F(2:end-1); % Solve the linear system
                                 % for interior nodes

end
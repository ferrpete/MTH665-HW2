% GStiff.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the 1D global stiffness matrix for FEM.

function A = GStiff(x)

% ===Input: vector x of mesh nodes===
% ===Output: Assembled stiffness matrix A===

N = length(x)-1; % number of elements
A = zeros(N+1,N+1); % initialize the stiffness matrix to zero

for i = 1:N % loop over elements
    h = x(i+1)-x(i); % element length
    n = [i i+1];
    A(n,n) = A(n,n) + [1 -1; -1 1]/h; % incorporate local stiffness
                                       % matrix into global matrix
end

end
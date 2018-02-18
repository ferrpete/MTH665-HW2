function A = GStiffRef(x)
% ===Input: vector x of mesh nodes ===
% ===Output: Assembled stiffness matrix A ===

N = length(x) - 1; % number of elements
A = zeros(N+1, N+1); % initialize the stiffness matrix to zero

for i=1:N % loop over elements
    
    xi = [x(i),x(i+1)]; % element nodes
    Aloc = StiffALoc(xi); % local stiffness matrix
    n = [i i+1]; % local to global map
    A(n,n) = A(n,n) + Aloc; % Local to Global
    
end
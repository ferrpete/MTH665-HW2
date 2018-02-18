function A = GStiffRef(x)
% ===Input: vector x of mesh nodes ===
% ===Output: Assembled stiffness matrix A ===

N = (length(x) - 1)/2; % number of elements
A = zeros(2*N+1, 2*N+1); % initialize the stiffness matrix to zero

for i=1:N % loop over elements
    
    xi = [x(2*i-1),x(2*i),x(2*i+1)]; % element nodes
    Aloc = StiffALoc(xi); % local stiffness matrix
    n = [2*i-1 2*i 2*i+1]; % local to global map
    A(n,n) = A(n,n) + Aloc; % Local to Global
    
end
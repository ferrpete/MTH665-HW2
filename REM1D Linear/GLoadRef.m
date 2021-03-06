% GLoadRef.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the 1D load vector for FEM.

function F = GLoadRef(x)

% ===Input: vector x of mesh nodes===
% ===Output: Load vector F using Trapezoidal Rule===

N = length(x)-1;
F = zeros(N+1,1);

for i = 1:N
    h = x(i+1) - x(i);
    n = [i, i+1];
    F(n) = F(n) + LoadLoc([x(i),x(i+1)]);
    
end

end
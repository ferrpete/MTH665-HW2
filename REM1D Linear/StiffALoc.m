function Aloc = StiffALoc(xi)
% ===Input xi = [x(i) x(i+1)], local nodes for element ===
% ===Output Aloc = Local stiffness matrix for element ===

Aloc = zeros(2,2); % initialize local matrix
jac = (xi(2) - xi(1))/2; % jacobian of map

[qWts, qPts] = Trapezoidal(-1,1);

for i = 1:length(qWts) % loop over quadrature points
    
    x = xi + (1 + qPts(1))*jac; % x on physical element = to qPts on ref.
    [psi, dpsi] = Reference(qPts(1)); % evaluate shape function
    Aloc = Aloc + dpsi/jac*dpsi'/jac*qWts(1)*jac;
    
end
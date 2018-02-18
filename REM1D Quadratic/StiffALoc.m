function Aloc = StiffALoc(xi)
% ===Input xi = [x(i) x(i+1)], local nodes for element ===
% ===Output Aloc = Local stiffness matrix for element ===

Aloc = zeros(3,3); % initialize local matrix
jac = (xi(3) - xi(1))/2; % jacobian of map

[qWts, qPts] = Simpson(-1,1);

for i = 1:length(qWts) % loop over quadrature points
    
    x = xi + (1 + qPts(1))*jac; % x on physical element = to qPts on ref.
    [psi, dpsi] = Reference(qPts(i)); % evaluate shape function
    Aloc = Aloc + dpsi/jac*dpsi'/jac*qWts(i)*jac;
    
end
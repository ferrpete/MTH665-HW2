function Floc = LoadLoc(xi)
% ===Input xi = [x(i) x(i+1)], local nodes for element ===
% ===Output Floc = Local load vector for element ===

Floc  = zeros(3,1); % initialize local load vector
jac = (xi(3) - xi(1))/2; % jacobian of map

[qWts,qPts] = Simpson(-1,1);

for i=1:length(qWts) % loop over quadrature points
    
    x = xi + (1+qPts(1))*jac; % x on physical element = to qPts on ref.
    [psi, dpsi] = Reference(qPts(i)); % evaluate shape function
    fx = feval(@Loadf,x); % Load is a user defined function for f.
    Floc = Floc + fx'.*psi*qWts(i)*jac;
    
end
function [psi, dpsi] = Reference(x)

psi(1,1) = 0.5 - 0.5*x;
psi(2,1) = 0.5 + 0.5*x;

dpsi(1,1) = -0.5;
dpsi(2,1) = 0.5;

end
function [psi, dpsi] = Reference(x)

psi(1,1) = (x.*(x+1))./2;
psi(2,1) = -((x+1).*(x-1));
psi(3,1) = (x.*(x-1))./2;

dpsi(1,1) = (2*x + 1)./2;
dpsi(2,1) = -2*x;
dpsi(3,1) = (2*x - 1)./2;

end
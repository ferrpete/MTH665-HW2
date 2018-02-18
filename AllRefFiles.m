% MainRef.m
% Peter Ferrero, Oregon State University, MTH655, 2/17/2018
% The main file for the 1D reference element method to solve Problem 4 of Homework 2 for
% MTH 655.

clear all

n = [2,4,8];
N = length(n);
x = [0:0.001:1]';
ExactSol = Exact(x);

figure(1)
plot(x, ExactSol, 'k')
xlabel('x')
ylabel('u(x)')
hold on

for i = 1:N
    
    [FemSol, x] = SimpleFEM1DRef(n(i));
    ExactSol = Exact(x');
    plot(x,FemSol,'-o')
    
    dExactSol = pi.*cos(pi.*x');
    errorMax(i) = norm(ExactSol-FemSol, inf);
    errorL2(i) = sqrt(Simpson13Approx(n(i),x,(ExactSol-FemSol).^2));
    fL2(i) = sqrt(Simpson13Approx(n(i),x,Loadf(x).^2));
    
    for j = 2:length(FemSol)
        
        dError(j-1) = ((FemSol(j) - FemSol(j-1))/(x(j) - x(j-1))) - dExactSol(j-1);
        
    end
    
    errorE(i) = sqrt(Simpson13Approx(n(i),x(2:end),dError.^2));
    
end

a = 0; % left endpoint
b = 1; % right endpoint
h = (b-a)./n; % uniform mesh size

legend('Exact', 'h = 0.5', 'h = 0.25', 'h = 0.125')
hold off

figure(2)
loglog(h,h,'k--',h,h.^2,'k-',h,errorMax,'*-r',h,errorL2,'*-b',h,errorE,'*-c')
xlabel('Mesh size, h', 'Interpreter', 'latex')
ylabel('Error, e(h)', 'Interpreter', 'latex')
legend({'Linear', 'Quadratic', '$\max_i |u(x_i)-u_k(x_i)|$', '$\left \| u - u_k \right \|_{L^2}$', '$\left \| u - u_k \right \|_E$'}, 'Interpreter', 'latex')
legend('Location', 'southeast')

% SimpleFEM1DRef.m
% Peter Ferrero, Oregon State University, MTH655, 1/31/2018
% A simple 1D reference element method to solve Problem 4 of Homework 2 for MTH 655.

function [FemSol, x] = SimpleFEM1DRef(N)

a = 0; % left endpoint
b = 1; % right endpoint
h = (b-a)/N; % uniform mesh size
x = a:h:b; % mesh nodes

A = GStiffRef(x); % Global Stiffness matrix
F = GLoadRef(x); % Load vector

FemSol = zeros(N+1, 1); % Initialize the FEM solution
FemSol(2:N) = A(2:N,2:N)\F(2:N); % Solve the linear system
                                 % for interior nodes

end

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

end

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

end

function [psi, dpsi] = Reference(x)

psi(1,1) = 0.5 - 0.5*x;
psi(2,1) = 0.5 + 0.5*x;

dpsi(1,1) = -0.5;
dpsi(2,1) = 0.5;

end

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

function Floc = LoadLoc(xi)
% ===Input xi = [x(i) x(i+1)], local nodes for element ===
% ===Output Floc = Local load vector for element ===

Floc  = zeros(2,1); % initialize local load vector
jac = (xi(2) - xi(1))/2; % jacobian of map

[qWts,qPts] = Trapezoidal(-1,1);

for i=1:length(qWts) % loop over quadrature points
    
    x = xi + (1+qPts(1))*jac; % x on physical element = to qPts on ref.
    [psi, dpsi] = Reference(qPts(i)); % evaluate shape function
    fx = feval(@Loadf,x); % Load is a user defined function for f.
    Floc = Floc + fx'.*psi*qWts(i)*jac;
    
end

end

% Loadf.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the local 1D load vector for FEM.

function f = Loadf(x)

% ===Input = x, mesh nodes at which to evaluate the load function f
% ===Output = f, value of load f at x

f = (pi^2)*sin(pi*x);

end
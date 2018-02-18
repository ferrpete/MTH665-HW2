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
% Main.m
% Peter Ferrero, Oregon State University, MTH655, 2/13/2018
% The main file for the FEM 1D method to solve Problem 1 of Homework 2 for
% MTH 655.

clear all

n = [2,4,8,16,32,64,128,256,512,1024];
N = length(n);
x = [0:0.001:1]';
ExactSol = Exact(x);

figure(1)
plot(x, ExactSol, 'k')
xlabel('x')
ylabel('u(x)')
hold on

for i = 1:N
    
    [FemSol, x] = SimpleFEM1D(n(i));
    ExactSol = Exact(x');
    errorMax(i) = norm(ExactSol-FemSol, inf);
    errorL2(i) = sqrt(Simpson13Approx(n(i),x,(ExactSol-FemSol).^2));
    fL2(i) = sqrt(Simpson13Approx(n(i),x,Loadf(x).^2));
    dError(1) = 0;
    
    for j = 2:length(FemSol)
        
        dError(j) = ((ExactSol(j) - FemSol(j)) -...
            (ExactSol(j-1) - FemSol(j-1)))/(x(j) - x(j-1));
        
    end
    
    errorE(i) = sqrt(Simpson13Approx(n(i),x,dError.^2));
    plot(x,FemSol,'-o')
    
end

a = 0; % left endpoint
b = 1; % right endpoint
h = (b-a)./n; % uniform mesh size

legend('Exact', 'h = 0.5', 'h = 0.25', 'h = 0.125')
hold off

figure(2)
loglog(h,h.^2,'k-',h,errorMax,'*-r',h,errorL2,'*-b',h,errorE,'*-c')
xlabel('Mesh size, h', 'Interpreter', 'latex')
ylabel('Error, e(h)', 'Interpreter', 'latex')
legend({'Quadratic', '$\max_i |u(x_i)-u_k(x_i)|$', '$\left \| u - u_k \right \|_{L^2}$', '$\left \| u - u_k \right \|_E$'}, 'Interpreter', 'latex')
legend('Location', 'southeast')

figure(3)
loglog(h,errorE,'*-c',h,h.*fL2,'k-')
title('Convergence of FEM Solution in Energy Norm', 'Interpreter', 'latex')
xlabel('Mesh size, h','Interpreter','latex')
legend({'$\left \| u - u_k \right \|_E$', '$h \left \| f \right \|_{L^2}$'}, 'Interpreter', 'latex')
legend('Location', 'southeast')

figure(4)
loglog(h,errorL2,'*-b',h,(h.^2).*fL2,'k-')
title('Convergence of FEM Solution in $L_2$ Norm', 'Interpreter', 'latex')
xlabel('Mesh size, h','Interpreter','latex')
legend({'$\left \| u - u_k \right \|_{L^2}$', '$h^2 \left \| f \right \|_{L^2}$'}, 'Interpreter', 'latex')
legend('Location', 'southeast')
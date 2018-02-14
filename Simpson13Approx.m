function [approx] = Simpson13Approx(n, x, Integrand)

h = (x(end)-x(1))/n; % Interval length
approx = (h/3)*(Integrand(1) + Integrand(end));

for i=2:n
    
    if rem(i,2) == 0
    
        approx = approx + 4*(h/3)*Integrand(i);
        
    else
        
        approx = approx + 2*(h/3)*Integrand(i);
        
    end
    
end

end
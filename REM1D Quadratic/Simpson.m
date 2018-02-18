function [qWts, qPts] = Simpson(c,d)
% qPts are the endpoints of the reference element

qPts(1) = c;
qPts(2) = (c+d)/2;
qPts(3) = d;

qWts(1) = (d-c)/6;
qWts(2) = 4*((d-c)/6);
qWts(3) = (d-c)/6;

end
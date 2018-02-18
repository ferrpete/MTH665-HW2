function [qWts, qPts] = Trapezoidal(c,d)
% qPts are the endpoints of the reference element

qPts(1) = c;
qPts(2) = d;

qWts(1) = (d-c)/2;
qWts(2) = (d-c)/2;

end
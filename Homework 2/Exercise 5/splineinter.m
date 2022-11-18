function [coef] = splineinter(x0,y0,dy0,x1,y1,dy1)

% obtain polynomial coefficients a, b, c, d
% a + b(x-x0) + c(x-x0)^2 + d(x-x0)^3

a = y0;
b = dy0;
dx =x1-x0;
cd = [3/dx^2 -1/dx; -2/dx^3 1/dx^2]*[y1-y0-dy0*dx; dy1-dy0];
c = cd(1);
d = cd(2);
coef = [a b c d];
clear all;
%close all hidden;
clc;

N = 10;
aint = 2;
yint = [0 2 4 6];
alpha = zeros(N,1); alpha(1) = 1;
beta  = zeros(N,1); beta(1)  = 1;
gamma = zeros(N,1); gamma(1) = 1;

[thetaest,X,C,V] = viterbiangleestimation(aint,yint,alpha,beta,gamma);
disp(thetaest);


draw_diagram(X,C,V);
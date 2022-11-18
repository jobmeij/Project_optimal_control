clear all;
close all hidden;
clc;

N = 10;
aint = 3;
yint = [0 2 6 8];
alpha = zeros(N,1); alpha(1) = 1/2; alpha(2) = 1/4; alpha(N) = 1/4;
beta  = zeros(N,1); beta(1)  = 1/2; beta(2)  = 1/4; beta(N)  = 1/4;
gamma = zeros(N,1); gamma(1) = 1/2; gamma(2) = 1/4; gamma(N) = 1/4;

[thetaest, X, V, C] = viterbiangleestimation_job(aint,yint,alpha,beta,gamma);
disp(thetaest);


draw_diagram_job(X,C,V);
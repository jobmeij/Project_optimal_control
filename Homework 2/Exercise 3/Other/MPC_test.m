clear all, close all, clc

optpolicy = 4;
fdisturbances = 0;
% double integrator model
Ac = [0 1;0 0];
Bc = [0;1];
tau = 0.2;
Q = [1 0; 0 1];
S = [0;0];
R = 0.01;
QT = [1 0; 0 1];
n = 2;
m = 1;
x0 = [1;0];
c = 0.25;
sigmaw = 0.01; % disturbance level
H = 25; % prediction horizon for mpc
h = 50; % simulation horizon
sysd = c2d(ss(Ac,Bc,zeros(1,n),0),tau);
A = sysd.a; B = sysd.b; C = sysd.c;

%% Preliminaries
% preliminaries for the policies
switch optpolicy
case 1
P{h+1} = QT;
for k = h:-1:1 % Riccati equations
P{k} = A'*P{k+1}*A + Q - (S+A'*P{k+1}*B)*pinv(R+B'*P{k+1}*B)*(S'+B'*P{k+1}*A);
K{k} = -pinv(R+B'*P{k+1}*B)*(S'+B'*P{k+1}*A);
end
case 2
U2 = quadconstrainedcontrol(A,B,Q,R,c,h,x0,1);
case 4
U3 = [quadconstrainedcontrol(A,B,Q,R,c,H,x0,2);
zeros(h-H,1)];
end
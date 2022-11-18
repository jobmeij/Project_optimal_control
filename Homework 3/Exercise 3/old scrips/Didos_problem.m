% Dido’s problem
% 4SC000 optimal control and dynamic programming
% Homework 3, Exercise 3
% Marcel van Wensveen
% 22-01-2019

clear all, close all, clc

%
syms t x(t)u_x(t) lambda1 lambda2 gamma
u = [u_x; u_y];
lambda = [lambda1; lambda2];
u_x(t) = diff(x(t))

f = (-gamma*u./(sqrt(1+u.^2)) + lambda);
s = solve(f==0, lambda);
pretty(s.lambda1)
pretty(s.lambda2)

diff(
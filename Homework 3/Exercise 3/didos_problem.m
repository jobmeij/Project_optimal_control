clear all;
close all hidden;
clc;

%% Calculate system (y_dot)
syms y y_dot Lambda LT L0 x
syms yx(x) Lambdax(x)
Dyx = diff(yx);

A = [0];
B = [1];
f = A*y + B*y_dot; % y_dot

cost = y;
constraint = sqrt(1+y_dot^2);

Lagran = cost-Lambda*constraint;

Lambda_dot = -A';
Euler_lagrange_eq = y_dot*jacobian(Lagran,y_dot)-Lagran - c1;

sol = solve(Euler_lagrange_eq,y_dot);
u = sol(1); %y_dot
disp('y_dot(x) = ')
pretty(u)
%% Calculate constants and simulate system
clc
syms c1_sym c2_sym lambda_sym
y0 = 0;
yT = 0;
x0 = -2;
xT = 30;
Ropelength = 35;
dt = 0.001;
if (Ropelength < (xT-x0))
    disp('Rope lenghth is <= than the distance between the begin and end point, no solutions')
elseif (Ropelength > ((xT-x0)*pi/2))
    disp('Optimal curve with this rope length is larger than half a circle, unable to solve')
end

eq_begin = ((x0 + c2_sym)^2+(y0+c1_sym)^2 == lambda_sym^2);
eq_end = ((xT + c2_sym)^2+(yT+c1_sym)^2 == lambda_sym^2);
sol = solve([eq_begin; eq_end], [c2_sym, c1_sym]);
c2 = double(sol.c2_sym(end))
c2_syms_eq = sol.c1_sym(2)==c1_sym;
sol = solve(c2_syms_eq,lambda_sym);
Lambda2_eq = sol(1)^2;
sol = vpasolve(2*sqrt(Lambda2_eq)*atan(sqrt((x0+c2)^2)/c1_sym)==Ropelength,c1_sym);
c1 = double(sol(end))

lambda_solution = sqrt((x0 + c2)^2 + (y0+c1)^2);
disp(['The value of gamma is: ',num2str(lambda_solution)])
x_vec = x0:dt:xT;
y_vec = zeros(1,length(x_vec));
ropelength_vec = zeros(1,length(x_vec));
y_dot_vec = zeros(1,length(x_vec));
for i = 1:length(x_vec)
    y_vec(i) = sqrt((x0+c2)^2+(y0+c1)^2-(x_vec(i)+c2)^2)-c1;
end   
for i=2:length(y_vec)
    y_dot_vec(i-1) = (y_vec(i)-y_vec(i-1))/dt;
    ropelength_vec(i) = ropelength_vec(i-1)+ sqrt(1+ y_dot_vec(i-1)^2)*dt;
end

figure(1)
plot(x_vec,y_vec)
title('Rope in state space')
xlabel('x')
ylabel('y(x)')
grid on
ylim([0, (xT-x0)/2])
figure(2)
plot(x_vec, ropelength_vec)
title('Cummulative sum of rope length')
xlabel('x')
ylabel('Length')
grid on
figure(3)
plot(x_vec, y_dot_vec)
title('ydot')
xlabel('x')
grid on
ylabel('ydot(x)')

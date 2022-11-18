
clear all;
close all hidden;
clc;

Length = 3.5;
syms L1 L2 G

Eq1 = ((3)/sqrt(G^2-L1^2))-3;
Eq2 = ((3)/sqrt(G^2-L2^2));
Eq3 = 3*sqrt(1+((L2^2)/(G^2-L2^2))^2) - Length;

S = solve([Eq1 == 0 Eq2 == 0 Eq3 == 0],[L1 L2 G]);
S
%%

sol = 1;
T = 3;
dt = 0.001;
t = 0:dt:T;
L = zeros(2,length(t)+1);
L(:,1) = [double(S.L1(sol)); double(S.L2(sol))];
gamma = double(S.G(sol));
x = zeros(2,length(t)+1);
u = zeros(2,length(t));
x(:,1) = [0;0];

for i = 1:length(t)
    u(:,i) = [sqrt(L(1,i)^2/(gamma^2-L(1,i)^2));sqrt(L(2,i)^2/(gamma^2-L(2,i)^2))];
    L(:,i+1) = L(:,i);% + [0;1]*dt;
    x(:,i+1) = x(:,i) + u(:,i)*dt;
end
figure(1)
plot(t,L(1,1:end-1),t,L(2,1:end-1))
title('lambda')
legend('\lambda_1','\lambda_2')
figure(2)
plot(x(1,:),x(2,:));
title('state space x')
xlabel('x')
ylabel('y')
figure(3)
plot(t,u(1,:),t,u(2,:))
title('u')
legend('u_x','u_y')
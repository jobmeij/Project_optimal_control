%% 4SC000 Homework 2 Exercise 1
clear all; close all; clc;

% Run learner solution.
Tend  = 2;
N = 3;
n = 3;
T = [0:N-1]/(N-1)*Tend;
rho = 0.01;
P = [0 0 0;
    1 1 0;
    4 3 0]';
tau1 = 0.5;
X0 = [P(:,1) (P(:,2)-P(:,1))/T(2) zeros(3,1)];
N12 = 10;
[xix,xiy,xiz]  = lqrtrajgeneration(P,T,X0,rho,tau1,N12)

% Run reference solution.
[xix_,xiy_,xiz_]  = lqrtrajgeneration(P,T,X0,rho,tau1,N12)

plot3(P(1,:),P(2,:),P(3,:))
hold on
plot3(xix(1,:),xiy(1,:),xiz(1,:))
axis equal
grid on
xlabel('x')
ylabel('y')
zlabel('z')
view(-90,90)
hold off

ex = xix(:)-xix_(:);
ey = xiy(:)-xiy_(:);
ez = xiz(:)-xiz_(:);

yReference = 1;
if norm(ex,inf) < 0.1 & norm(ey,inf) < 0.1 & norm(ez,inf) < 0.1
    y = 1;
else
    y = 0;
end



% Compare.
assessVariableEqual('y', yReference);
clear all;
close all hidden;
clc;

% Run learner solution.
Tend  = 2;
N = 3;
n = 3;
T = [0:N-1]/(N-1)*Tend;
rho = 0.01;
P = [0 0 0;1 1 0;4 3 0]';
tau1 = 0.5;
X0 = [P(:,1) (P(:,2)-P(:,1))/T(2) zeros(3,1)];
N12 = 10;

[xix,xiy,xiz]  = lqrtrajgeneration(P,T,X0,rho,tau1,N12);

figure();
plot3(P(1,:),P(2,:),P(3,:),".-k");
hold on;
plot3(xix(1,:),xiy(1,:),xiz(1,:),".b");
plot3(xix(1,1),xiy(1,1),xiz(1,1),".g");
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
view(0,90);
hold off;

% figure(2);
% plot3(xix(2,:),xiy(2,:),xiz(2,:),".-k");
% hold on;
% plot3(xix(2,1),xiy(2,1),xiz(2,1),".-g");
% grid on;
% view(0,90);
% xlabel('x');
% ylabel('y');
% zlabel('z');
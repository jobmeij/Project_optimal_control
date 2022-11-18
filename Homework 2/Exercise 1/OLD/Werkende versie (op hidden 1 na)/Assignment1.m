clear all;
close all hidden;
clc;

thetavec = (0:1:4)/5*(2*pi);
thetavec2 = thetavec+1/(10)*(2*pi);
p = [];
for i = 1:5
   p = [p;
   2*sin(thetavec(i)) 2*cos(thetavec(i));
   sin(thetavec2(i))  cos(thetavec2(i))];
end
p = [p;p(1,:)];
P = [p';zeros(1,size(p,1))];
Tend  = 20;
N = size(P,2);
n = 3;
T = [0:N-1]/(N-1)*Tend;
rho = 0.01;
tau1 = 0.1;
X0 = [P(:,1) (P(:,2)-P(:,1))/T(2) zeros(3,1)];%
N12 = 10;

[xix,xiy,xiz]  = lqrtrajgeneration(P,T,X0,rho,tau1,N12);

figure();
plot3(P(1,:),P(2,:),P(3,:),".-k");
hold on;
plot3(xix(1,:),xiy(1,:),xiz(1,:),"-b");
plot3(xix(1,1),xiy(1,1),xiz(1,1),"og");
plot3(xix(1,end),xiy(1,end),xiz(1,end),"or");
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
view(0,90);

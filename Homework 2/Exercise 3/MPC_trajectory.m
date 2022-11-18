%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimal control and dynamic programming    %
%Homework 2                                 %
%Assignment 3                               %
%Date: 09-Jan-2019                          %
%Group: 2                                   %
%Bob Clephas            | 1271431           %
%Tom van de laar        | 1265938           %
%Job Meijer             | 1268155           %
%Marcel van Wensveen    | 1253085           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all hidden;
clc;

rho = 0.01;
tau = 0.1;
H = 80;
c = 0.4;
N = 400;
x0 = [10 -20 -1;
      0  -10 -1;
      0  0   0.1;
      0  0   0.4];
t_vect = (0:tau:N*tau)';
x = zeros(4,length(t_vect));
u = zeros(1,length(t_vect));
r = 3*sin(0.3*t_vect);

A = [0 1 0 0;
     0 0 1 0
     0 0 0 1
     0 0 0 0];
B = [0;0;0;1];
C = [1 0 0 0];
D = 0;

sys = ss(A,B,C,D);
sysd = c2d(sys,tau);
Ad = sysd.A;
Bd = sysd.B;
Cd = C;

PHImpc = zeros(H,4);
for i = 1:H
    PHImpc(i,:) = Cd*Ad^(i); 
end
RHOmpc = zeros(H,H);
for i = 1:H
    h = Cd*Ad^(i-1)*Bd;
    for j = 1:H-i+1
        RHOmpc(j+i-1,j) = h;
    end
end

Q = 1;
Qmpc = eye(H)*Q;

R = rho;
Umpc = eye(H)*R;

Hmpc = RHOmpc'*Qmpc*RHOmpc + Umpc;

LB = -ones(H,1)*c;
UB = ones(H,1)*c;
%A_bound_mpc = [eye(H);-eye(H)];
%B_bound_mpc = [c*ones(H,1);c*ones(H,1)];
options = optimoptions('quadprog','Display','off');

for init = 1:3
x(:,1) = x0(:,init);
for t = 1:N+1
    Rmpc = 3*sin(0.3*(t+1:t+H)'*tau);
    bmpc = Rmpc - PHImpc*x(:,t);
    Gmpc = -RHOmpc'*Qmpc*bmpc;
    [U,J] = quadprog(Hmpc,Gmpc,[],[],[],[],LB,UB,[],options);
%    [U,J] = quadprog(Hmpc,Gmpc,A_bound_mpc,B_bound_mpc,[],[],[],[],[],options);
    u(t) = U(1);
    x(:,t+1) = Ad*x(:,t) + Bd*u(t);
end

fig = figure(init+2);
subplot(311);
plot(t_vect,r,t_vect,x(1,1:(N+1)));
plotname = ['Init conditions set ', num2str(init)];
fig.Name = (plotname);
fig.NumberTitle = 'off';
title('Position P_k')
grid on;
xlabel("Time [k\tau] [s]");
ylabel("Position");
legend('Reference','Position P_k')
subplot(312);
plot(t_vect,x(2,1:(N+1)));
title('Velocity v_k')
grid on;
xlabel("Time [k\tau] [s]");
ylabel("Velocity");
subplot(313);
plot(t_vect,u);
title('Control input u_k')
grid on;
xlabel("Time [k\tau] [s]");
ylabel("Input");
ylim([-c c]);
end
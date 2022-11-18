%Demonstration of the Robin Sharp Optimal Preview Steering Control
%Alexander Brown
%March 15, 2012
clear all
close all
clc

%Adapted from Sharp, valetsiotis, 2001


%vehicle parameters
m = 1573;
Iz = 2550;
b = 1.73;
a = .913;
Cf = 2*88310;
Cr = 2*64076;
G = 16;%ratio of hand wheel to road wheel angle.
g = 9.8;


U = 38.9; %forward speed, m/s
T = .02;%sampling time

%vehicle state vector is [y ydot psi psidot]'

%Set up the continuous-time vehicle dynamics
A = [0 1 0 0;...
    0 -(Cf+Cr)/(m*U) (Cf+Cr)/m (b*Cr-a*Cf)/(m*U);...
    0 0 0 1;...
    0 (b*Cr-a*Cf)/(Iz*U) (a*Cf-b*Cr)/Iz -(a^2*Cf+b^2*Cr)/(Iz*U)];

B = [0 Cf/(G*m) 0 a*Cf/(G*Iz)]';%^sharp uses that G term here to scale, but i see no need at thsi point.

sysd = c2d(ss(A,B,eye(4),[0 0 0 0]'),T); %discrete time state space description.

%break
%this set a b c d variables.

%this is in ERROR coordinates

%now we decide how many preview points to use
np = 250;%try ten preview points first.
D = zeros(np);%initialize
D(1:end-1,2:end) = eye(np-1);%set up the shift register
E = zeros(np,1);%input vector for shift register
E(end) = 1;%

%now we need to set up the cost function
C = zeros(2,np+4);%initialize to the right size
C(:,1:6) = [1 0 0 0 -1 0;0 0 1 0 1/(U*T) -1/(U*T)];%set actual values
Q = [.1 0; 0 0];%weighting on heading vs lateral position error
R1 = C'*Q*C;
R2 = 1;

%this is our total cost function for our entire system. We need to only
%look at the system without any preview. 

Ad = sysd.a;
Bd = sysd.b;

Aaug = blkdiag(Ad,D);%automatically creates the big a matrix
Baug = [Bd;zeros(size(E))];%input vector for steering angle
B2aug = [zeros(size(Bd));E];%input vector for road

%now we have to set up the road.

% Uroad = zeros(1,sim_steps);
% Uroad(250:end) = 3.5;
% Uroad(750:end) = 0;
% Uroad(1:249) = 0;%set up some initial error in the road

Uroad = [zeros(1,2/T) linspace(0,3,4/T) 3*ones(1,10/T)];
%lag the reference vector
Uroad_lagged = [Uroad(1)*ones(1,np) Uroad(1:end-np)]
sim_steps = length(Uroad);
t = T*(1:sim_steps)-T;

X0 = zeros(size(Baug));%set initial conditions to zero

%try finding stuff for full aug sys
[Kf,Sf,Ef] = dlqr(Aaug,[Baug],C'*Q*C,[1]);%steering only

%try finding a normal LQR just acting on position error
C2 = [1 0 0 0 -1; 0 0 0 0 0];
Ad2 = [Ad [0;0;0;0];[0 0 0 0 0]];
Bd2 = [Bd;0];
[Klqr,Sflqr,Eflqr] = dlqr(Ad2,Bd2,C2'*Q*C2,[1]);

%different C matrix for simulation. Define output as y only
C_sim = Kf;
C_sim2 = zeros(1,length(B2aug));
C_sim2(1)=1;
Clqr2 = [1 0 0 0 0];
Clqr = Klqr;


%now let's see what happens when we simulate the system
[y,x] = dlsim(Aaug-[Baug]*Kf,[B2aug],C_sim,0,Uroad,X0);
[ylqr,xlqr] = dlsim(Ad2-Bd2*Klqr,[.5/18.5;0;0;0;0],Clqr,0,Uroad_lagged,[0;0;0;0;0]);%simulate normal lqr

sys_lqr = ss(Ad2-Bd2*Klqr,[.5/18.5;0;0;0;0],Clqr2,0,T);
sys = ss(Aaug-[Baug]*Kf,[B2aug],C_sim2,0,T);

figure
bode(sys,sys_lqr);
legend('Preview Controller','LQR');

title(['Preview LQR vs LQR: Bandwidth ' num2str(bandwidth(sys)/(2*pi)) ' Hz vs ' num2str(bandwidth(sys_lqr)/(2*pi)) ' Hz']);

% 
% figure
% bode(sys_lqr)
% title(['bode plot for non-preview LQR: Bandwidth ' num2str(bandwidth(sys_lqr)/(2*pi)) ' Hz'])%now, Uroad is actually the offset at the FARTHEST POINT AHEAD THAT WE
% %MEASURE!!



%now let's plot car lateral position
figure
plot(t,Uroad_lagged,'r',t,x(:,1)','k',t,xlqr(:,1),'g')
title('lateral position as a function of time')
legend('road centerline','car lateral position (preview)','car lateral position (LQR)')
xlabel('time (s)')
ylabel('global y-coordinate (m)')

%now let's plot car steering angle position
figure
plot(t,-y','k',t,ylqr,'g')
title('steering angle as a function of time')
legend('steering angle (preview)','steering angle (LQR)')
xlabel('time (s)')
ylabel('roadwheel angle (rad)')

%now let's plot the preview gains
figure
plot(Kf(5:end),'k.')
xlabel('preview point number')
ylabel('preview gain')
title('preview gain vs. preview point number')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimal control and dynamic programming    %
%Homework 3                                 %
%Assignment 4                               %
%Date: 28-Jan-2019                          %
%Group: 2                                   %
%Bob Clephas            | 1271431           %
%Tom van de laar        | 1265938           %
%Job Meijer             | 1268155           %
%Marcel van Wensveen    | 1253085           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all hidden;
clc;

A = [0 1;-1 0]; % [0 1;-1 0];
B = [0; 4]; % [0; 4]
x0_vect = [8 -9 2;0 -8 -1];
x0_init = 3;
umin = -1;
umax = 1;
Maxswitching = 10; % 2
dt = 0.001;
Nmax = 30000;
Stepvect = [1 10 100];
tol = 0.05;

syms s;
Svect = 0:0.1:pi;
Sumin = zeros(2,length(Svect));
Sumax = zeros(2,length(Svect));

for u = [umin umax]
    X = int((expm(A*-s))*B*u,s);
    for i = 1:length(Svect)
        t = Svect(i);
        x = double(subs(X,s,0))-double(subs(X,s,t));
        switch u
            case umax
                Sumax(:,i) = x;
            case umin
                Sumin(:,i) = x;
        end
    end
end

x0 = x0_vect(:,x0_init);

bestdist = inf;
bestT = inf;
bestX = [];
besta = 0;
Switching = 0;
Switchingdone = false;
bestSP = [];


figure(2);
set(2,'Position',[900 556 500 430]);
figure(3);
set(3,'Position',[1400 556 500 430]);
figure(4);
set(4,'Position',[400 556 500 430]);

while ~Switchingdone
    for i = Stepvect
        a_vect = besta-(pi/i):0.1/i:besta+(pi/i);
        for a = a_vect
            L1 = cos(a);
            L2 = sin(a);
            L = zeros(2,Nmax);
            L(:,1) = [L1; L2];
            x = zeros(2,Nmax+1);
            x(:,1) = x0;
            u = zeros(1,Nmax);
            H = zeros(1,Nmax);
            SP = [];

            count = 0;
            done = false;
            Ht = inf;
            curbestdist = inf;
            curbestT = 0;
            t = 1;
            while t < Nmax && ~done && (L(1,t) ~= 0 || L(2,t) ~= 0)
                u(t) = -umax*sign(L(:,t)'*B);
                if u(t) == 0
                    u(t) = umax;
                end
                if t > 1 && u(t) ~= u(t-1)
                    count = count + 1;
                    SP = [SP t];
                end
                
                xdot = A*x(:,t) + B*u(t);
                x(:,t+1) = x(:,t) + xdot*dt;
                
                H(t) = L(:,t)'*(A*x(:,t)+B*u(t)) + 1;
                prevHt = Ht;
                Ht = H(t);
                Ldot = -A'*L(:,t);
                L(:,t+1) = L(:,t) + Ldot*dt;
                
                if count == Switching+1
                    done = true;
                end
                if  norm(x(:,t+1)) <= curbestdist
                    curbestdist = norm(x(:,t+1));
                    curbestT = t;
                end
                t = t + 1;
            end
            t = curbestT;
            if norm(x(:,t+1)) <= bestdist
                bestdist = norm(x(:,t+1));
                bestL = L(:,1:t);
                bestH = H(:,1:t);
                bestX = x(:,1:t+1);
                bestU = u(:,1:t);
                bestN = t;
                bestT = t*dt;
                besta = a;
                bestSP = SP(SP<t);
                
                
                figure(2);
                hold off;
                plot(bestX(1,:),bestX(2,:),'Linewidth',2);
                hold on;
                plot(bestX(1,1),bestX(2,1),'.g','markerSize',20);
                plot(bestX(1,t+1),bestX(2,t+1),'.r','markerSize',20);
                plot(Sumin(1,:),Sumin(2,:),'.r','markerSize',10);
                plot(Sumax(1,:),Sumax(2,:),'.g','markerSize',10);
                for i = 1:size(bestSP,2)
                    plot(bestX(1,bestSP(i)),bestX(2,bestSP(i)),'.c','markerSize',20);
                end
                
                bv = 1:2:(Switching+1)*2;
                for b = bv
                    plot(4*b-4*cos(Svect), -4*sin(Svect),'.g','markerSize',10);
                    plot(-4*b-4*cos(Svect), 4*sin(Svect),'.r','markerSize',10);
                end
                
                plot(0,0,'k.','markerSize',25);
                grid minor;
                title("State space L0 = ["+L1 +";"+L2+"] N = "+bestN);
                
                figure(3);
                subplot(211);
                plot((0:t-1)*dt,bestH,'.','Markersize',10);
                grid minor;
                title("Hamiltonian");
                subplot(212);
                plot((0:t-1)*dt,abs(bestL'*B),'.','Markersize',10);
                grid minor;
                title("|L'*B|");
                
                figure(4);
                plot(bestT,bestdist,'.b','Markersize',20);
                xlabel("T");
                ylabel("Dist");
                grid on;
                hold on;
                title("Convergence");
                drawnow();
            end
        end
    end
    if ~isempty(bestX) && norm(bestX(:,end)) < tol || Switching == Maxswitching
        Switchingdone = true;
    else
        Switching = Switching + 1;
    end
end
disp("Done");

Xmin = min(bestX(1,:))-5;
Xmax = max(bestX(1,:))+5;


figure(1);
set(1,'Position',[980 70 900 900]);
subplot(312);
plot((0:bestN)*dt,bestX(1,:),'Linewidth',2);
grid minor;
title("Position");
xlabel('time t')
subplot(313);
plot((0:bestN)*dt,bestX(2,:),'Linewidth',2);
grid minor
title("Velocity");
xlabel('time t')
subplot(311);
plot((0:bestN-1)*dt,bestU,'Linewidth',2);
grid minor;
title("Input T = "+bestT);
xlabel('time t')

figure(2);
set(2,'Position',[50 70 900 900]);
hold off;
plot(bestX(1,:),bestX(2,:),'Linewidth',2);
hold on;

bv = 1:2:(Switching+1)*2;
Svect = 0:0.05:pi;
for b = bv
    plot(4*b-4*cos(Svect), -4*sin(Svect),'.g','markerSize',5);
    plot(-4*b-4*cos(Svect), 4*sin(Svect),'.r','markerSize',5);
end

plot(bestX(1,1),bestX(2,1),'.g','markerSize',20);
plot(bestX(1,bestN+1),bestX(2,bestN+1),'.r','markerSize',20);

for i = 1:size(bestSP,2)
    plot(bestX(1,bestSP(i)),bestX(2,bestSP(i)),'.c','markerSize',20);
end

title('State space solution');
legend('Optimal solution','Switching plane u = 1;','Switching plane u = -1');
xlabel('Position');
ylabel('Velocity');

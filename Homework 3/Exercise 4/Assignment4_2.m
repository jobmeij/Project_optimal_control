clear all;
close all hidden;
clc;

A = [0 1;-1 0]; % [0 1;-1 0];
B = [0; 4]; % [0; 4]
x0_vect = [8 -29 2;0 -8 -1];
x0_init = 2;
umin = -1;
umax = 1;
dt = 0.001;
Nmax = 5000;
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
x(:,1) = x0;
u = zeros(1,Nmax);
L = zeros(2,Nmax);
L(:,1) = [-1;-1];
SP = [];
for t = 1:Nmax
    u(t) = -umax*sign(L(:,t)'*B);
    if u(t) == 0
        u(t) = umax;
    end
    u(t) = umax;
    if t > 1 && u(t) ~= u(t-1)
        SP = [SP t];
    end
    xdot = A*x(:,t) + B*u(t);
    x(:,t+1) = x(:,t) + xdot*dt;
    Ldot = -A'*L(:,t);
    L(:,t+1) = L(:,t) + Ldot*dt;
end

figure(1);
set(1,'Position',[1400 42 500 430]);
subplot(311);
plot((0:Nmax)*dt,x(1,:),'Linewidth',2);
grid minor;
title("Position T = "+Nmax*dt);
subplot(312);
plot((0:Nmax)*dt,x(2,:),'Linewidth',2);
grid minor;
title("Velocity");
subplot(313);
plot((0:Nmax-1)*dt,u,'Linewidth',2);
grid minor;
title("Input");

figure(2);
set(2,'Position',[900 556 500 430]);
plot(x(1,:),x(2,:),'Linewidth',2);
hold on;
plot(x(1,1),x(2,1),'.g','markerSize',20);
plot(x(1,Nmax),x(2,Nmax),'.r','markerSize',20);
plot(Sumin(1,:),Sumin(2,:),'.r','markerSize',10);
plot(Sumax(1,:),Sumax(2,:),'.g','markerSize',10);
for i = 1:size(SP,2)
    plot(x(1,SP(i)),x(2,SP(i)),'.c','markerSize',20);
end
%plot(4*sin(Svect), 4-4*cos(Svect),'.g','markerSize',10);
%plot(-4*sin(Svect), -4+4*cos(Svect),'.r','markerSize',10);
plot(0,0,'k.','markerSize',25);
grid minor;
title("State space");

figure(3);
set(3,'Position',[1400 556 500 430]);
plot((0:t)*dt,abs(L'*B),'.','Markersize',10);
grid minor;
title("|L'*B|");


Xmin = min(x(1,:))-5;
Xmax = max(x(1,:))+5;
figure(5);
set(5,'Position',[900 42 500 430]);
fps = 30;
for t = 1:ceil((1/fps)/dt):size(x,2)
    figure(5);
    hold off;
    plot(x(1,t),0,'sb','Markersize',50);
    hold on;
    plot([x(1,t) x(1,t)+u(t)],[0 0],'k','linewidth',3);
    plot([x(1,t) x(1,t)+x(2,t)],[0 0],'r','linewidth',2);
    plot(0,0,'.r','Markersize',20);
    plot(x0(1),0,'.g','Markersize',20);
    grid on;
    title("Animation");
    xlim([Xmin Xmax])
    ylim([-10 10]);
    drawnow();
    pause(1/fps);
end





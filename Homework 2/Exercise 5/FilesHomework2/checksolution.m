clear all;
close all hidden;
clc;

% draw race track ---------------------------------------------------------
figure(1);
set(1,'Position',[0  40  1920  960]);
subplot(1,2,1);
drawracecircuit;
% -------------------------------------------------------------------------

indx0 = 2; % 1 | 2
optinput = 1; %1 | 2
save_u = false;
switch optinput
    case 1
        if indx0 == 1
            StartPos = [942 822]';
        else
            StartPos = [717 964]';
        end
        plot(StartPos(1),StartPos(2),'ow');
        Reference = createReference(Track,30);
        [u] = mintimeAV(Reference, Track, StartPos, Finish);
        if save_u
            save(strcat("Generated_u",string(indx0),".mat"),'u');
        end
    case 2
        if indx0 == 1
            load Generated_u1;
        else
            load Generated_u2;
        end
end
if indx0 == 1
    opponend = load('u1');
else
    opponend = load('u2');
end
optanimation = 1; % 1 makes animation, 0 does not make animation
% -------------------------------------------------------------------------

% define model-------------------------------------------------------------
tau  = 0.1;
Ac   = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
Bc   = [0 0 1 0;0 0 0 1]';
Cc   = [1 0 0 0;0 1 0 0];
Dc   = [0 0;0 0];
sysd = c2d(ss(Ac,Bc,Cc,Dc),tau);
A = sysd.a; B = sysd.b;
n = size(A,2);
% -------------------------------------------------------------------------

% Calculate opponend ------------------------------------------------------
h_ = size(opponend.u,2);
x_ = zeros(n,h_+1);
if indx0 == 1
    x0 = [942 822]';
else
    x0 = [717 964]';
end
x_(1:2,1) = x0;
for k=1:h_
    x_(:,k+1) = A*x_(:,k)+B*opponend.u(:,k);
end
% -------------------------------------------------------------------------

% Calculate model ---------------------------------------------------------
h = size(u,2);
x = zeros(n,h+1);
if indx0 == 1
    x0 = [942 822]';
else
    x0 = [717 964]';
end
x(1:2,1) = x0;
for k=1:h
    x(:,k+1) = A*x(:,k)+B*u(:,k);
end
% -------------------------------------------------------------------------

% plot and animations -----------------------------------------------------
h = min(h,h_);
figure(1);
hold on;
plot(x0(1),x0(2),'o','color',[1 0 0]);
xf = x(1:2,k+1);
line([896 977],[622 622],'LineWidth',2,'color',[0.5 0.5 0.5]) % finish line
axis([0 1024 0 1100]);
xlabel('x'); ylabel('y');
subplot(1,2,1);
if optanimation == 1
    H1_1 = patch([-10 -10 10 10],[-10 10 10 -10],[1 1 1]);
    H1_2 = patch([-10 -10 10 10],[-10 10 10 -10],[1 0 0]);
    t1 = hgtransform('Parent',gca);
    t2 = hgtransform('Parent',gca);
    set(H1_1,'parent',t1);
    set(H1_2,'parent',t2);
    for k=1:h+1
        T_1 = [ [eye(3) [x_(1,k) x_(2,k) 0]'];zeros(1,3) 1];
        T_2 = [ [eye(3) [x(1,k) x(2,k) 0]'];zeros(1,3) 1];
        set(t1,'Matrix',T_1);
        set(t2,'Matrix',T_2);
        pause(tau);
    end
end
for i = 1:size(x,2)
    switch OnTrack(Track,x(:,i))
        case true
            figure(1);
            plot(x(1,i),x(2,i),'.g','Markersize',15);
        case false
            figure(1);
            plot(x(1,i),x(2,i),'.r','Markersize',15);
    end
end
plot(xf(1),xf(2),'dw');
title(sprintf('Trajectory time: %0.2f',tau*h));
subplot(1,2,2);
plot(tau*[0:h-1],sqrt(u(1,:).^2+u(2,:).^2));
xlabel('time');
ylabel('norm control input');
% -------------------------------------------------------------------------
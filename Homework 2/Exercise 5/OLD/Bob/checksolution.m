clear all;
close all hidden;
clc;

% draw race track ---------------------------------------------------------
figure(1);
set(1,'Position',[300  300  950  400]);
subplot(1,2,1);
drawracecircuit;
% -------------------------------------------------------------------------

indx0 = 0;
%  change to one and provide your function mintimeAV ----------------------
optinput = 1;
if(optinput == 0) 
    if indx0 == 1
        load u1;
    else
        load u2;
    end
else
    if indx0 == 1
        StartPos = [942 822]';
    else
        StartPos = [717 964]';
    end
    plot(StartPos(1),StartPos(2),'ow');
    Reference = createreference(Track, StartPos);
    [u] = mintimeAV(Reference, Track, StartPos, Finish);
end
optanimation = 0; % 1 makes animation, 0 does not make animationx
% -------------------------------------------------------------------------


% define model-------------------------------------------------------------
tau  = 0.1;
Ac   = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
Bc   = [0 0 1 0;
        0 0 0 1]';
Cc =   [1 0 0 0;
        0 1 0 0];
Dc = [0 0;
      0 0];
sysd = c2d(ss(Ac,Bc,Cc,Dc),tau);
A = sysd.a; B = sysd.b;
n = size(A,2);
% -------------------------------------------------------------------------

% run model-------------------------------------------------------------
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
    t1 = hgtransform('Parent',gca);
    set(H1_1,'parent',t1);
    for k=1:h+1
        T_2 = [ [eye(3) [x(1,k) x(2,k) 0]'];zeros(1,3) 1];
        set(t1,'Matrix',T_2);
        pause(tau);
    end
end
for i = 1:size(x,2)
    switch OnTrack(Track,x(:,i))
        case true
            figure(1);
            plot(x(1,i),x(2,i),'.g','Markersize',20);
        case false
            figure(1);
            plot(x(1,i),x(2,i),'.r','Markersize',20);
    end
end
plot(xf(1),xf(2),'dw');
title(sprintf('Trajectory time: %0.2f',tau*h));
subplot(1,2,2);
plot(tau*[0:h-1],sqrt(u(1,:).^2+u(2,:).^2));
xlabel('time');
ylabel('norm control input');
% -------------------------------------------------------------------------
close all, clear all, clc

indx0 = 1;
%  change to one and provide your function mintimeAV ----------------------
optinput = 0;
if(optinput == 0) 
    if indx0 == 1
        load u1;
    else
        load u2;
    end
else
    [u] = mintimeAV;
end
optanimation = 1; % 1 makes animation, 0 does not make animation
% -------------------------------------------------------------------------

% draw race track ---------------------------------------------------------
figure(1)
subplot(1,2,1)
drawracecircuit
% -------------------------------------------------------------------------

% define model-------------------------------------------------------------
tau  = 0.1;
Ac   = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
Bc   = [0 0 1 0;
        0 0 0 1]';
sysd = c2d(ss(Ac,Bc,[1 0 0 0],0),tau);
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
hold on
plot(x0(1),x0(2),'o','color',[1 0 0]);
xf = x(1:2,k+1);
line([896 977],[622 622],'LineWidth',2,'color',[0.5 0.5 0.5]) % finish line
axis([0 1024 0 1100])
xlabel('x'); ylabel('y')
subplot(1,2,1)
if optanimation == 1
    H1_1 = patch([-10 -10 10 10],[-10 10 10 -10],[1 1 1])
    t1 = hgtransform('Parent',gca);
    set(H1_1,'parent',t1)
    for k=1:h+1
        T_2 = [ [eye(3) [x(1,k) x(2,k) 0]'];zeros(1,3) 1];
        set(t1,'Matrix',T_2)
        pause(tau)
    end
end
plot(x(1,:),x(2,:),'color',[ 0 1 0])
plot(xf(1),xf(2),'d','color',[1 1 1]);
title('red circle- initial condition, white diamond final position')
text(100,20,sprintf('Trajectory time: %0.2f',tau*h),'color',[1 1 1]);
subplot(1,2,2)
plot(tau*[0:h-1],sqrt(u(1,:).^2+u(2,:).^2))
set(1,'Position',[97.4  207  900  420])
xlabel('time')
ylabel('norm control input')
% -------------------------------------------------------------------------
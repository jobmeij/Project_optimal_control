function u = mintimeAV(Reference, Track, Startpos, Finish);
u = zeros(2,1);
% define model-------------------------------------------------------------
tau  = 0.1;
Ac   = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
Bc   = [0 0 1 0;0 0 0 1]';
Cc   = [1 0 0 0;0 1 0 0];
Dc   = [0 0;0 0];
sysd = c2d(ss(Ac,Bc,Cc,Dc),tau);
Ad = sysd.a;
Bd = sysd.b;
Cd = sysd.c;
% -------------------------------------------------------------------------
F = 1000;
rho = 0;
Umax = 300;
H = 30;
Hbound = 5;
Nlaps = 1;
Laps = -1;
N = size(Track,2);
Nref = size(Reference,2);

PHImpc = zeros(H*2,4);
for i = 1:H
    PHImpc(i*2-1:i*2,:) = Cd*Ad^(i);
end

RHOmpc = zeros(H*2,H*2);
for i = 1:H
    h = Cd*Ad^(i-1)*Bd; %2x2
    for j = 1:H-i+1
        RHOmpc((j+i)*2-3:(j+i)*2-2,j*2-1:j*2) = h;
    end
end

Q = 1;
Qmpc = eye(H*2)*Q;

R = rho;
Umpc = eye(H*2)*R;

Hmpc = RHOmpc'*Qmpc*RHOmpc + Umpc;
Hmpc = (Hmpc+Hmpc')/2; %make matrix symmetric

resolution = 16;
n = resolution/4;
Ainputmpc = zeros(resolution*H,2*H);
for i = 1:H %Create input bounds
    for a = 1:n
        angle0 = (pi/2)/n*(a-1);
        angle1 = (pi/2)/n*a;
        x0 = cos(angle0);
        y0 = sin(angle0);
        x1 = cos(angle1);
        y1 = sin(angle1);
        dx = abs(x1-x0);
        dy = abs(y1-y0);
        Ainputmpc(a+(i-1)*resolution,i*2-1:i*2) = [dx dy];
        Ainputmpc(a+n+(i-1)*resolution,i*2-1:i*2) = [dx -dy];
        Ainputmpc(a+n*2+(i-1)*resolution,i*2-1:i*2) = [-dx dy];
        Ainputmpc(a+n*3+(i-1)*resolution,i*2-1:i*2) = [-dx -dy];
    end
end
Binputmpc = ones(resolution*H,1)*Umax*Ainputmpc(n,1);

Aoutputmpc0 = zeros(Hbound*2,H*2);
Aoutputmpc1 = zeros(Hbound*2,H*2);
for i = 1:Hbound
    Aoutputmpc0(i,:) = RHOmpc(i*2-1,:)+RHOmpc(i*2,:);
    Aoutputmpc0(i+Hbound,:) = RHOmpc(i*2-1,:)+RHOmpc(i*2,:);
end
Boutputmpc0 = zeros(Hbound*2,1);
Boutputmpc1 = zeros(Hbound*2,1);

Ub = ones(H*2,1)*(Umax);
Lb = ones(H*2,1)*(-Umax);
options = optimoptions('quadprog','Display','off');

% run model-------------------------------------------------------------
Rmpc = zeros(H*2,1);
x = zeros(4,2);
x(1:2,1) = Startpos;
k = 0;
prevpart = 1;
inside = true;
figure(2);
set(2,'Position',[960 40 960 960]);
figure(3);
set(3,'Position',[0 40 960 960]);
U = zeros(H*2,1);

while Laps < Nlaps && k < F && inside
    k = k+1;
    part = trackpart(Track,x(1:2,k),prevpart);
    prevpart = part;
    prevparth = part;
    figure(2);
    hold off;
    plot([Track(1,:) Track(1,1)],[Track(2,:) Track(2,1)],'.-b');
    hold on;
    plot([Track(3,:) Track(3,1)],[Track(4,:) Track(4,1)],'.-b');
    plot(x(1,k),x(2,k),'.g','Markersize',20);
    
    trackindex = 0;
    bestdist = inf;
    for r = 1:Nref
        dist = norm(Reference(1:2,r)-x(1:2,k));
        if dist < bestdist
            trackindex = r;
            bestdist = dist;
        end
    end
    
    for i = 1:H
        trackindex = mod(trackindex,Nref)+1;
        Rmpc(i*2-1:i*2) = Reference(1:2,trackindex);
        plot(Rmpc(i*2-1),Rmpc(i*2),'.r','Markersize',10);
    end
    
    for h = 1:Hbound
        xest = PHImpc(h*2-1,:)*x(:,k) + RHOmpc(h*2-1,:)*U;
        yest = PHImpc(h*2,:)*x(:,k)   + RHOmpc(h*2,:)*U;
        plot(xest,yest,'.m','Markersize',20);
        cartpos = [xest yest];
        part = trackpart(Track,cartpos,prevparth);
        prevparth = part;
        
        dxtrack(1) = Track(1,mod(part,N)+1)-Track(1,part);
        dytrack(1) = Track(2,mod(part,N)+1)-Track(2,part);
        dxtrack(2) = Track(3,mod(part,N)+1)-Track(3,part);
        dytrack(2) = Track(4,mod(part,N)+1)-Track(4,part);

        for i = 1:2
            dx = dxtrack(i);
            dy = dytrack(i);
            if(abs(dx) >= abs(dy))
                dydx = dy/dx;
                norm1 = sqrt(1+dydx^2);
                yoffset = (Track(i*2,part)-(Track(i*2-1,part)*dydx));
                cartoffset = PHImpc(h*2,:)*x(:,k)-(PHImpc(h*2-1,:)*x(:,k)*dydx);
                if(yoffset+dydx*Track(5,part)>=Track(6,part))
                    Aoutputmpc1(h,:) = repmat([-dydx 1]./norm1,1,H);
                    Boutputmpc1(h) = yoffset/norm1;
                    Boutputmpc0(h) = cartoffset/norm1;
                else
                    Aoutputmpc1(Hbound+h,:) = -repmat([-dydx 1]./norm1,1,H);
                    Boutputmpc1(Hbound+h) = -yoffset/norm1;
                    Boutputmpc0(Hbound+h) = -cartoffset/norm1;
                end
            else
                dxdy = dx/dy;
                norm1 = sqrt(1+dxdy^2);
                xoffset = (Track(i*2-1,part)-(Track(i*2,part)*dxdy));
                cartoffset = PHImpc(h*2-1,:)*x(:,k)-(PHImpc(h*2,:)*x(:,k)*dxdy);
                if(xoffset+dxdy*Track(6,part)>=Track(5,part))
                    Aoutputmpc1(h,:) = repmat([1 -dxdy]./norm1,1,H);
                    Boutputmpc1(h) = xoffset/norm1;
                    Boutputmpc0(h) = cartoffset/norm1;
                else
                    Aoutputmpc1(Hbound+h,:) = -repmat([1 -dxdy]./norm1,1,H);
                    Boutputmpc1(Hbound+h) = -xoffset/norm1;
                    Boutputmpc0(Hbound+h) = -cartoffset/norm1;
                end
            end
        end
    end
    %drawbounds(Aoutputmpc1,Boutputmpc1,2,Hbound,0,1024,0,1100,['r','k']);
    xlabel("x");
    ylabel("y");
    title("MPC");
    axis([0 1024 0 1100]);
    
    Aoutputmpc = Aoutputmpc0 .* Aoutputmpc1;
    Boutputmpc = -(Boutputmpc0 -  Boutputmpc1);
    
    Ampc = [Ainputmpc; Aoutputmpc];
    Bmpc = [Binputmpc; Boutputmpc];
    Fmpc = -RHOmpc'*Qmpc*(Rmpc - PHImpc*x(:,k));
    
    U = quadprog(Hmpc,Fmpc,Ampc,Bmpc,[],[],Lb,Ub,[],options);
    
    if size(U,1) == 2*H
        for h = 1:H
            xest = PHImpc(h*2-1,:)*x(:,k) + RHOmpc(h*2-1,:)*U;
            yest = PHImpc(h*2,:)*x(:,k)   + RHOmpc(h*2,:)*U;
            plot(xest,yest,'.c','Markersize',15);
        end
        u(:,k) = U(1:2);
        figure(3);
        plot(u(1,k),u(2,k),'ro');
        hold on;
        s = 0:0.01:1;
        Ubound = Umax*exp(2*pi*s*sqrt(-1));
        plot(real(Ubound),imag(Ubound),'g');
        %drawbounds(Ainputmpc,Binputmpc,resolution,1,-Umax,Umax,-Umax,Umax,['b']);
        xlabel("x");
        ylabel("y");
        title("Bounded model input");
        axis([-Umax Umax -Umax Umax]);
        x(:,k+1) = Ad*x(:,k)+Bd*u(:,k);
        if passfinish(Finish,x(1:2,k),x(1:2,k+1))
            Laps = Laps + 1;
        end
        inside = OnTrack(Track,x(1:2,k+1));
    else
        inside = false;
        disp("No solution found!");
    end
    drawnow();
end
end
function u = mintimeAV(Reference, Track, Startpos, Finish);
u = zeros(2,1);
% define model-------------------------------------------------------------
tau  = 0.1;
Ac   = [0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0];
Bc   = [0 0 1 0;
        0 0 0 1]';
Cc   = [1 0 0 0;
        0 1 0 0];
Dc   =  [0 0;
         0 0];
sysd = c2d(ss(Ac,Bc,Cc,Dc),tau);
Ad = sysd.a; 
Bd = sysd.b;
Cd = sysd.c;
% -------------------------------------------------------------------------
F = 200;
rho = 0.5;
Umax = 300;
H = 20;
Nlaps = 1;
Laps = -1;
N = size(Track,2);

PHImpc = zeros(H*2,4);
for i = 1:H
    PHImpc(i*2-1:i*2,:) = Cd*Ad^(i-1);
    %PHImpc(i*2-1:i*2,:) = Cd*Ad^(i); 
end

RHOmpc = zeros(H*2,H*2);
for i = 1:H
    h = Cd*Ad^(i-1)*Bd; %2x2
    %h = Cd*Ad^(i)*Bd; % 2x2
    for j = 1:H-i+1
        RHOmpc((j+i)*2-3:(j+i)*2-2,j*2-1:j*2) = h;
    end
end

Q = 1;
Qmpc = eye(H*2)*Q;

R = rho;
Umpc = eye(H*2)*R;

Hmpc = RHOmpc'*Qmpc*RHOmpc + Umpc;

Ainputmpc = zeros(4*H,2*H);
for i = 1:H
    Ainputmpc(1+(i-1)*4,i*2-1:i*2) = [1 1];
    Ainputmpc(2+(i-1)*4,i*2-1:i*2) = [1 -1];
    Ainputmpc(3+(i-1)*4,i*2-1:i*2) = [-1 1];
    Ainputmpc(4+(i-1)*4,i*2-1:i*2) = [-1 -1];
end

Binputmpc = ones(4*H,1)*Umax;

Aoutputmpc0 = zeros(2,H*2);
Aoutputmpc0(1,:) = RHOmpc(1,:)+RHOmpc(2,:);
Aoutputmpc0(2,:) = -(RHOmpc(1,:)+RHOmpc(2,:));

Boutputmpc = zeros(2,1);

options = optimoptions('quadprog','Display','off');

% run model-------------------------------------------------------------

Rmpc = zeros(H*2,1);
x = zeros(4,2);
x(1:2,1) = Startpos;
k = 0;
prevpart = 1;
inside = true;
Npoints = 1;
prevpoint = [];
figure(2);
set(2,'Position',[1000  100  400  400]);
%figure(3);
%set(3,'Position',[1450  100  400  400]);
while Laps < Nlaps && k < F && inside
    k = k+1;
    part = trackpart(Track,x(1:2,k),prevpart);
    prevpart = part;
    
    figure(2);
    hold off;
    plot([Track(1,:) Track(1,1)],[Track(2,:) Track(2,1)],'.-b');
    hold on;
    plot([Track(3,:) Track(3,1)],[Track(4,:) Track(4,1)],'.-b');
    plot(x(1,k),x(2,k),'.g','Markersize',10);
    maxdist = 2;
    trackindex = part;
    if k == 1
        prevpoint = x(1:2,k);
    end
%     for i = 1:H
%         if Npoints >= 1
%             prevpoint = Track(5:6,trackindex);
%             trackindex = mod(trackindex + 1,N)+1;
%         end
%         diff = Track(5:6,trackindex)-prevpoint;
%         trackdist = sqrt(diff(1)^2+diff(2)^2)
%         Npoints = ceil(trackdist / maxdist)
%         refpoint = prevpoint+diff./Npoints;
%         prevpoint = refpoint;
%         Rmpc(i*2-1:i*2) = refpoint;
%         plot(Rmpc(i*2-1),Rmpc(i*2),'.r','Markersize',10);
%     end
    for i = 1:H
        trackindex = mod(trackindex,N)+1;
        Rmpc(i*2-1:i*2) = Track(5:6,trackindex);
        plot(Rmpc(i*2-1),Rmpc(i*2),'.r','Markersize',10);
    end
    dxoutside = (Track(1,part)-Track(1,mod(part+1,N)+1));
    dyoutside = (Track(2,part)-Track(2,mod(part+1,N)+1));
    dxinside = (Track(3,part)-Track(3,mod(part+1,N)+1));
    dyinside = (Track(4,part)-Track(4,mod(part+1,N)+1));
    
    for i = 1:2
        switch i
            case 1
                dx = dxoutside;
                dy = dyoutside;
            case 2
                dx = dxinside;
                dy = dyinside;
        end
        if(abs(dx) >= abs(dy))
            dydx = dy/dx;
            yoffset = Track(i*2,part)-(Track(i*2-1,part)*dydx);
            line([0 1024],[yoffset yoffset+1024*dydx]);
        else
            dxdy = dx/dy;
            xoffset = Track(i*2-1,part)-(Track(i*2,part)*dxdy);
            line([xoffset xoffset+1100*dxdy],[0 1100]);
        end
    end
    axis([0 1024 0 1100]);
    %drawnow();
    
    for i = 1:2
        switch i
            case 1
                dx = dxoutside;
                dy = dyoutside;
            case 2
                dx = dxinside;
                dy = dyinside;
        end
        if(abs(dx) >= abs(dy))
            %disp("dydx");
            dydx = dy/dx;
            angle = atan(dydx);
            yoffset = Track(i*2,part)-(Track(i*2-1,part)*dydx);
            if(yoffset+dydx*x(1,k)>=x(2,k))
                Aoutputmpc(1,:) = Aoutputmpc0(1,:).*repmat([dy dx],1,H);
                Boutputmpc(1) = yoffset-PHImpc(2,:)*x(:,k);
            else
                Aoutputmpc(2,:) = Aoutputmpc0(2,:).*repmat([dy dx],1,H);
                Boutputmpc(2) = -(yoffset-PHImpc(2,:)*x(:,k));
            end
        else
            %disp("dxdy");
            dxdy = dx/dy;
            angle = atan(dxdy);
            xoffset = Track(i*2-1,part)-(Track(i*2,part)*dxdy);
            if(xoffset+dxdy*x(2,k)>=x(1,k))
                Aoutputmpc(1,:) = Aoutputmpc0(1,:).*repmat([dy dx],1,H);
                Boutputmpc(1) = xoffset-PHImpc(1,:)*x(:,k);
            else
                Aoutputmpc(2,:) = Aoutputmpc0(2,:).*repmat([dy dx],1,H);
                Boutputmpc(2) = -(xoffset-PHImpc(1,:)*x(:,k));
            end
        end
    end
    %Aoutputmpc
    %Boutputmpc
    %figure(3);
    %hold off;
    if(Aoutputmpc(1,1) == Aoutputmpc0(1,1))
        rc = Aoutputmpc(1,2)/Aoutputmpc(1,1);
        plot([Boutputmpc(1) 1024*rc+Boutputmpc(1)],[0 1100],'r');
    else
        rc = Aoutputmpc(1,1)/Aoutputmpc(1,2);
        plot([0 1024],[Boutputmpc(1) 1100*rc+Boutputmpc(1)],'r');
    end
    hold on;
    if(Aoutputmpc(2,1) == Aoutputmpc0(2,1))
        rc = Aoutputmpc(2,2)/Aoutputmpc(2,1);
        plot([Boutputmpc(2) 1024*rc+Boutputmpc(2)],[0 1100],'b');
    else
        rc = Aoutputmpc(2,1)/Aoutputmpc(2,2);
        plot([0 1024],[Boutputmpc(2) 1100*rc+Boutputmpc(2)],'b');
    end
    axis([0 1024 0 1100]);
    drawnow();
    
    %Ampc = [Ainputmpc; Aoutputmpc];
    %Bmpc = [Binputmpc; Boutputmpc];
    Ampc = [Ainputmpc];
    Bmpc = [Binputmpc];
    Fmpc = -RHOmpc'*Qmpc*(Rmpc - PHImpc*x(:,k));
    
    U = quadprog(Hmpc,Fmpc,Ampc,Bmpc,[],[],[],[],[],options);
    if size(U,1) > 2
        u(:,k) = U(1:2);
        x(:,k+1) = Ad*x(:,k)+Bd*u(:,k);
        if passfinish(Finish,x(1:2,k),x(1:2,k+1))
            Laps = Laps + 1;
        end
        inside = OnTrack(Track,x(1:2,k+1));
        inside = true;
    else
        inside = false;
    end
    pause(tau);
end

end
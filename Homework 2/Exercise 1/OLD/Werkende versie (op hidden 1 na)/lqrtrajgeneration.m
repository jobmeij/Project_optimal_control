function [xix, xiy, xiz] = lqrtrajgeneration(P,T,X0,rho,tau1,N12)

N = length(T);%Number of input points
h = floor(T(end)/tau1);%Number of samples "discrete" path
n = size(X0,2);%Number of derivatives

tau2 = tau1/N12;%"continuous" path sampling period
H = N12*h+1;%Number of samples "continuous" path

r = zeros(3,h+1);%Sampled input path
u = zeros(3,h);%nth derivative

xix = zeros(n+1,H);
xiy = zeros(n+1,H);
xiz = zeros(n+1,H);

for i = 1:N-1%interpolate input path
    for t = T(i):tau1:T(i+1)
        k = round((t/tau1)+1);
        r(:,k) = P(:,i) + (((t - T(i))/(T(i+1)-T(i))) * (P(:,i+1)-P(:,i)));
    end
end

%State space model
Ax = zeros(n,n);
Ax(1:n-1,2:n) = eye(n-1,n-1);
Ay = Ax;
Az = Ax;
A = [Ax zeros(n,n) zeros(n,n);zeros(n,n) Ay zeros(n,n);zeros(n,n) zeros(n,n) Az];

Bx = zeros(n,1);
Bx(n) = 1;
By = Bx;
Bz = Bx;
B = [Bx zeros(n,1) zeros(n,1);zeros(n,1) By zeros(n,1);zeros(n,1) zeros(n,1) Bz];

Cx = zeros(1,n);
Cx(1) = 1;
Cy = Cx;
Cz = Cx;
C = [Cx zeros(1,n) zeros(1,n);zeros(1,n) Cy zeros(1,n);zeros(1,n) zeros(1,n) Cz];

D = zeros(3,3);

sys = ss(A,B,C,D);

%Discretisation
sysd = c2d(sys,tau1);
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

hp = h+1;
Bd0 = [Bd;zeros(hp*3,3)];
Cd0 = [Cd zeros(3,hp*3)];

%Reference shift register
Arx = zeros(hp,hp);
Arx(1:hp-1,2:hp) = eye(hp-1,hp-1);
Ary = Arx;
Arz = Arx;
Ar = [Arx zeros(hp,hp) zeros(hp,hp);zeros(hp,hp) Ary zeros(hp,hp);zeros(hp,hp) zeros(hp,hp) Arz];
Adr = [Ad zeros(n*3,hp*3); zeros(hp*3,n*3) Ar];

Brx = zeros(hp,1);
Brx(hp) = 1;
Bry = Brx;
Brz = Brx;
Br = [Brx zeros(hp,1) zeros(hp,1);zeros(hp,1) Bry zeros(hp,1);zeros(hp,1) zeros(hp,1) Brz];
Br0 = [zeros(n*3,3);Br];

Crx = zeros(1,hp);
Crx(1) = 1;
Cry = Crx;
Crz = Crx;
Cr = [Crx zeros(1,hp) zeros(1,hp);zeros(1,hp) Cry zeros(1,hp);zeros(1,hp) zeros(1,hp) Crz];
Cdr = [Cd -Cr];

%LQR
Q = Cdr'*eye(3)*Cdr;
R = eye(3,3)*rho;
S = zeros(3*(n+hp),3);

UC = zeros(3,3*(n+hp),h);
JC = zeros(3*(n+hp),3*(n+hp),h+1);
JC(:,:,h+1) = Q;

for k = h:-1:1
    JCkp = squeeze(JC(:,:,k+1));
    UC(:,:,k) = -(R+Bd0'*JCkp*Bd0)^-1*(Bd0'*JCkp*Adr+S');
    JC(:,:,k) = Adr'*JCkp*Adr - (Adr'*JCkp*Bd0+S)*(R+Bd0'*JCkp*Bd0)^-1*(Bd0'*JCkp*Adr+S') + Q;
end

%Input calculation
xdr = zeros(3*(n+hp),h);
xdr(:,1) = [X0(1,:)';X0(2,:)';X0(3,:)';r(1,:)';r(2,:)';r(3,:)'];
for k = 1:h
    UCk = squeeze(UC(:,:,k));
    u(:,k) = UCk*xdr(:,k);
    xdr(:,k+1) = Adr*xdr(:,k) + Bd0*u(:,k);
end

%Output calculation
sysd2 = c2d(sys,tau2);
Ad2 = sysd2.A;
Bd2 = sysd2.B;
Cd2 = sysd2.C;
Dd2 = sysd2.D;

xd2 = zeros(3*n,H);
xd2(:,1) = [X0(1,:)';X0(2,:)';X0(3,:)'];

for t = 1:H-1
    k = ceil(t/N12);
    if(N12 < 1)
        k = k-1;
    end
    xd2(:,t+1) = Ad2*xd2(:,t) + Bd2*u(:,k);
    xix(n+1,t) = u(1,k);
    xiy(n+1,t) = u(2,k);
    xiz(n+1,t) = u(3,k);
end

xix(1:n,:) = xd2(0*n+1:1*n,:);
xiy(1:n,:) = xd2(1*n+1:2*n,:);
xiz(1:n,:) = xd2(2*n+1:3*n,:);
xix(n+1,H) = xix(n+1,H-1);
xiy(n+1,H) = xiy(n+1,H-1);
xiz(n+1,H) = xiz(n+1,H-1);

figure();
plot3(r(1,:),r(2,:),r(3,:),'.r');
grid on;
view(0,90);
end
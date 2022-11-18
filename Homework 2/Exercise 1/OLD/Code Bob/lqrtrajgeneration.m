function [xix, xiy, xiz] = lqrtrajgeneration(P,T,X0,rho,tau1,N12)

N = length(T);%Number of input points
h = floor(T(end)/tau1);%Number of samples "discrete" path
n = size(X0,2)+1;%Number of derivatives + trajectory

tau2 = tau1/N12;%"continuous" path sampling period
H = N12*h+1;%Number of samples "continuous" path

r = zeros(h+1,3);%Sampled input path
p = zeros(h+1,3);%Sampled output path
u = zeros(h,3);%nth derivative

for i = 1:N-1%interpolate input path
    for t = T(i):tau1:T(i+1)
        k = uint16((t/tau1)+1);
        r(k,:) = P(:,i) + (((t - T(i))/(T(i+1)-T(i))) * (P(:,i+1)-P(:,i)));
    end
end

%"continuous" output path with derivatives
xix = zeros(n,H);
xiy = zeros(n,H);
xiz = zeros(n,H);

xix(1:n-1,1) = X0(1,:);
xiy(1:n-1,1) = X0(2,:);
xiz(1:n-1,1) = X0(3,:);

JC = zeros(h,3,3);
UC = zeros(h,2,3);
%!!!!!!!!!!!!!!!! Pk+1 != pk + uk !!!!!!!!!!!!!!!!!!!!!
UC(h,:,:) = [(2*r(h+1,:))./(2*(rho+1)); ...
    ones(1,3).*((-2)./(2*(rho+1)))];

UCh1 = squeeze(UC(h,1,:))';
UCh2 = squeeze(UC(h,2,:))';

JC(h,:,:) = [(rho+1)*UCh1.^2 - 2*r(h+1,:).*UCh1 + r(h,:).^2 + r(h+1,:).^2; ...
    2*UCh1 - 2*r(h,:) - 2*r(h+1,:) + (rho+1)*2*UCh1.*UCh2 - r(h+1,:).*UCh2; ...
    2*UCh2 + 2 + (rho+1)*UCh2.^2];

for k = h-1:-1:1
    JCkp1 = squeeze(JC(k+1,1,:))';
    JCkp2 = squeeze(JC(k+1,2,:))';
    JCkp3 = squeeze(JC(k+1,3,:))';
    
    UC(k,:,:) = [(-JCkp2)./(2*(rho+JCkp3)); ...
        (-2*JCkp3)./(2*(rho+JCkp3))];
    
    UCk1 = squeeze(UC(k,1,:))'; 
    UCk2 = squeeze(UC(k,2,:))';
    
    JC(k,:,:) = [r(k,:).^2 + rho*UCk1.^2 + JCkp3.*UCk1.^2 + JCkp2.*UCk1 + JCkp1; ...
        (2*UCk1.*UCk2).*(JCkp3+rho) - 2*r(k,:) + 2*UCk1 + JCkp2.*(1+UCk2); ...
        1 + rho*UCk2.^2 + JCkp3.*(1+UCk2.^2) + 2*UCk2];
end

UC
JC

% TO BE DETERMINED VIA LQR
u = [1.6040 0.8020 0;
     1.7035 0.8518 0;
     0.4838 0.2419 0;
     0.0374 0.0187 0]
%%%%%%%%%%%%%%%%%%%%%%%%%%

p(1,:) = X0(:,1);

for k = 1:h
    UCk1 = squeeze(UC(k,1,:))'; 
    UCk2 = squeeze(UC(k,2,:))';
    
    u(k,:) = (1/factorial(n-2)).*(UCk1 + UCk2 .* p(k,:)).*(tau1^(n-2));
        
    xix(n,(k-1)*N12+1) = u(k,1);
    xiy(n,(k-1)*N12+1) = u(k,2);
    xiz(n,(k-1)*N12+1) = u(k,3);
    
    for t = (N12*(k-1))+1:1:N12*k
        for d = 1:n % Calculate trajectory and derivatives
            for hd = d:n % derivative terms [x = x0 + x'*dt + (1/2)*x''*dt^2 ...]
                xix(d,t+1) = xix(d,t+1) + (1/factorial(hd-d))*xix(hd,t)*tau2^(hd-d);
                xiy(d,t+1) = xiy(d,t+1) + (1/factorial(hd-d))*xiy(hd,t)*tau2^(hd-d);
                xiz(d,t+1) = xiz(d,t+1) + (1/factorial(hd-d))*xiz(hd,t)*tau2^(hd-d);
            end
        end
    end
    p(k+1,:) = [xix(1,k*N12+1) xiy(1,k*N12+1) xiz(1,k*N12+1)];    
end

u

figure();
plot3(r(:,1),r(:,2),r(:,3),'.-r');
hold on;
plot3(P(1,:),P(2,:),P(3,:),'ok');
plot3(p(:,1),p(:,2),p(:,3),'*g');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
view(0,90);

end
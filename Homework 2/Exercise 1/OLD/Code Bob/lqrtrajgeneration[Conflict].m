function [xix, xiy, xiz, r] = lqrtrajgeneration(P,T,X0,rho,tau1,N12)

N = length(T);%Number of input points
h = floor(T(end)/tau1);%Number of samples "discrete" path
n = size(X0,2)+1;%Number of derivatives + trajectory

tau2 = tau1/N12;%"continuous" path sampling period
H = N12*h+1;%Number of samples "continuous" path

r = zeros(3,h+1);%Sampled input path
p = zeros(3,h+1);%Sampled output path
u = zeros(3,h);%nth derivative 

for i = 1:N-1%interpolate input path
    for t = T(i):tau1:T(i+1)-tau1
        k = (t/tau1)+1;
        r(:,k) = P(:,i) + (((t - T(i))/(T(i+1)-T(i))) * (P(:,i+1)-P(:,i)));
    end
end

%"continuous" output path with derivatives
xix = zeros(n,H);
xiy = zeros(n,H);
xiz = zeros(n,H);

xix(1:n-1,1) = X0(1,:);
xiy(1:n-1,1) = X0(2,:);
xiz(1:n-1,1) = X0(3,:);

% TO BE DETERMINED VIA LQR%
% x4 = [1.6040 1.7035 0.4838 0.0374];
% y4 = [0.8020 0.8518 0.2419 0.0187];
% z4 = [0 0 0 0];
%

% h = 4
C = zeros(h,2,3);
for k = h:-1:1
   C(k,:,:) = [(2*r(:,k+1))/(2*rho+2) ones(3,1)*(-2/(2*rho+2))]';    
end

C;

p = r(:,1);
for k = 1:h
    p;
    x4(k) = C(k,1,1) + C(k,2,1) * p(1);
    y4(k) = C(k,1,1) + C(k,2,1) * p(1);
    z4(k) = C(k,1,1) + C(k,2,1) * p(1);
    p = [x4(k); y4(k); z4(k)];
end
%


for k = 1:h
    xix(n,(k-1)*N12+1) = x4(k);
    xiy(n,(k-1)*N12+1) = y4(k);
    xiz(n,(k-1)*N12+1) = z4(k);
    
    for t = (N12*(k-1))+1:1:N12*k
        for d = 1:n % Calculate trajectory and derivatives
            for hd = d:n % derivative terms [x = x0 + x'*dt + (1/2)*x''*dt^2 ...]
                xix(d,t+1) = xix(d,t+1) + (1/factorial(hd-d))*xix(hd,t)*tau2^(hd-d);
                xiy(d,t+1) = xiy(d,t+1) + (1/factorial(hd-d))*xiy(hd,t)*tau2^(hd-d);
                xiz(d,t+1) = xiz(d,t+1) + (1/factorial(hd-d))*xiz(hd,t)*tau2^(hd-d);
            end
        end
    end
    
    %DO SOMETHING WITH THIS COST 
    p(:,k) = [xix(1,(k-1)*N12+1) xiy(1,(k-1)*N12+1) xiz(1,(k-1)*N12+1)];
    u(:,k) = [xix(n,(k-1)*N12+1) xiy(n,(k-1)*N12+1) xiz(n,(k-1)*N12+1)];
    Cost = norm(p(:,k)-r(:,k),2)^2 + rho*norm(u(:,k),2)^2;
end

end



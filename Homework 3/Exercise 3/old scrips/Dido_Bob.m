clear all;
close all hidden;
clc;

T = 3;
Length = 3.5;

syms G t L10 L20 L1T L2T x1 x2 u1 u2 L1 L2

x = [x1;x2];
u = [u1;u2];
L = [L1;L2];

A = [0 0;0 0];
B = [1 0;0 1];

Q = [0 0;0 0];
R = [0 0;0 1];
C = [0 1];

f = A*x + B*u
g = -(C*x+G*sqrt(1+C*u.^2))

Ldot = - jacobian(f,x)'*L - jacobian(g,x)'

eq1 = 0 == -jacobian(f,u)'*L + jacobian(g,u)'

Su1 = solve(eq1(1),u1,'Real',true)
Su2 = solve(eq1(2),u2,'Real',true)

n = 2;
if ~isempty(Su1) && ~isempty(Su2)
    V = [Su1(n)/L1 0; 0 Su2(n)/L2];
elseif ~isempty(Su1)
    V = [Su1(n)/L1 0; 0 0];
elseif ~isempty(Su2)
    V = [0 0; 0 Su2(n)/L2];
else
    V = [0 0; 0 0];
end

%V = [1/sqrt(G^2-L1T^2) 0; 0 1/sqrt(G^2-L2T^2)];
H = [A V;Q -A']

eH = (expm(H*T));
eH = subs(eH,[L1 L2],[L1T L2T])

Eq1 = eH*([0;0;L10;L20]) == [3;0;L1T;L2T];
Eq1 = simplify(Eq1)
if ~isempty(Su2)
    Eq2 = G*3*sqrt(1+(Su2(n))^2) == G*Length;
else
    Eq2 = G == G;
end
S = solve([Eq1' Eq2],[L10 L20 G],'Real',true);

disp(S.L10);
disp(S.L20);
disp(S.G);

if ~isempty(S.L10)
    sol = 1;
    T = 3;
    dt = 0.001;
    t = 0:dt:T;
    L = zeros(2,length(t)+1);
    L(:,1) = [double(S.L10(sol)); double(S.L20(sol))];
    gamma = double(S.G(sol));
    x = zeros(2,length(t)+1);
    u = zeros(2,length(t));
    x(:,1) = [0;0];
    
    for i = 1:length(t)
        u(:,i) = [sqrt(L(1,i)^2/(gamma^2-L(1,i)^2));sqrt(L(2,i)^2/(gamma^2-L(2,i)^2))];
        L(:,i+1) = L(:,i) + [0;1]*dt;
        x(:,i+1) = x(:,i) + u(:,i)*dt;
    end
    figure(1)
    plot(t,L(1,1:end-1),t,L(2,1:end-1))
    title('lambda')
    figure(2)
    plot(x(1,:),x(2,:));
    title('state space x')
    figure(3)
    plot(t,u(1,:),t,u(2,:))
    title('u')
else
    disp("No solution");
end

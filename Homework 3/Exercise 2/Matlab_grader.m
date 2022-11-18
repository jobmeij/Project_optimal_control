clear all, close all, clc

[Q,R,W,V] = gainslqgmargins;
tfplant = tf([1 0.3],[1 4])*tf([1 4],[1 -5])*tf([1],[1 -8]);
[A,B,C,D] = tf2ss(tfplant.num{1},tfplant.den{1});
n = size(A,1);
lqgcontroller = lqg(ss(A,B,C,0),blkdiag(Q,R),blkdiag(W,V));
[num,den] = ss2tf(lqgcontroller.a,lqgcontroller.b,lqgcontroller.c,0);
lqgcontrollertf = tf(-num,den);
vecgains = [1/2+0.01:0.1:20];
vecphase = exp(j*pi/180*[-40:0.5:40]);
flag = 1;
nyquist(tfplant*lqgcontrollertf)

for i = 1:length(vecgains)
    % gain
    %eig(minreal(-vecgains(i)*tfplant*lqgcontrollertf/(1+vecgains(i)*tfplant*lqgcontrollertf)))
    if any(eig(minreal(-vecgains(i)*tfplant*lqgcontrollertf/(1+vecgains(i)*tfplant*lqgcontrollertf)))>0)
        flag = 0;
    end
end

% Compare.
%assessVariableEqual('flag', 1);
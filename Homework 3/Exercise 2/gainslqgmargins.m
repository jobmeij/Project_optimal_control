function [Q,R,W,V] = gainslqgmargins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimal control and dynamic programming    %
%Homework 3                                 %
%Assignment 2                               %
%Date: 28-Jan-2019                          %
%Group: 2                                   %
%Bob Clephas            | 1271431           %
%Tom van de laar        | 1265938           %
%Job Meijer             | 1268155           %
%Marcel van Wensveen    | 1253085           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tfplant = tf([1 0.3],[1 4])*tf([1 4],[1 -5])*tf([1],[1 -8]);
[~,~,C,~] = tf2ss(tfplant.num{1},tfplant.den{1});

R = 1e-14;
Q = C'*C;
W = eye(2); 
V = eye(2);

end
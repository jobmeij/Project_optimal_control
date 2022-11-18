%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Homework 2 - Exercise 2
%   
%   Calling the function 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
test = 2;
switch test
     case 0 % example from assignment Cost ~13.15
         
        Ac = [0 1 ; -1 -1];
        tau = 0.1;
        Bc = [0;1];
        Cc = [1 0];
        sysd = c2d(ss(Ac,Bc,Cc,zeros(size(Cc,1),size(Bc,2))),tau);
        A = sysd.a;
        B = sysd.b;
        V = 0.01;
        W = eye(2);
        Q = eye(2);
        R = 0.01;
        C = Cc;
        h = 1;
        disp(['Case 0, with h = ',num2str(h)])
        cost = lqgmultirate(A,B,C,Q,R,W,V,h);
     case 1 % Test 1 from grader | Cost ~13.1348 / 13.1032 !12.99
        % Run learner solution.
        Ac = [0 1 ; -1 -1];
        tau = 0.1;
        Bc = [0;1];
        Cc = eye(2);
        sysd = c2d(ss(Ac,Bc,Cc,zeros(size(Cc,1),size(Bc,2))),tau);
        A = sysd.a;
        B = sysd.b;
        V = 0.01*eye(2);
        W = eye(2);
        Q = eye(2);
        R = 0.01*eye(1);
        C = Cc;
        h = 1;
        disp('Case 1')
        cost = lqgmultirate(A,B,C,Q,R,W,V,h);
     case 2 % Cost ~ 16.01
        Ac = [0 1; -4 -2];
        tau = 0.1;
        Bc = [0;1];
        Cc = [1 0];
        sysd = c2d(ss(Ac,Bc,Cc,0),tau);
        A = sysd.a;
        B = sysd.b;
        V = 1;
        W = eye(2);
        Q = eye(2);
        R = 0.001*eye(1);
        C = Cc;
        h = 7;
        disp('Case 2')
        cost = lqgmultirate(A,B,C,Q,R,W,V,h);
     case 3
        Ac = [0 1 0; 0 0 1; 4 -4 3];
        tau = 0.1;
        Bc = [0;0;1];
        Cc = [1 0 0];
        sysd = c2d(ss(Ac,Bc,Cc,0),tau);
        A = sysd.a;
        B = sysd.b;
        N = 1;
        V = 1;
        W = [1 0.2 0;0.2 1 0;0 0 10];
        Q = eye(3);
        R = 0.1*eye(1);
        C = Cc;
        h = 1;
        disp('Case 3')
        cost = lqgmultirate(A,B,C,Q,R,W,V,h);
end




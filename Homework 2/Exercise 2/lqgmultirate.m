function [cost]  = lqgmultirate(A,B,C,Q,R,W,V,h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimal control and dynamic programming    %
%Homework 2                                 %
%Assignment 2                               %
%Date: 09-Jan-2019                          %
%Group: 2                                   %
%Bob Clephas            | 1271431           %
%Tom van de laar        | 1265938           %
%Job Meijer             | 1268155           %
%Marcel van Wensveen    | 1253085           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Horizon
T = 100000; 

% Initialize variables
n_states = size(A,1);
n_outputs = size(C,1);
phi = zeros(n_states,n_states,T);
L = zeros(n_states,n_outputs,T);
cost_trace = zeros(T,1);

% Calculate LQR solution
[P_bar,~,~] = dare(A,B,Q,R);

% Calculate required matrix for trace method
X_bar = A'*P_bar*B*inv(R+B'*P_bar*B)*B'*P_bar*A;

for i = 1:T-1
    % prediction step
    phi(:,:,i+1) = A*phi(:,:,i)*A' + W;
    
    % Correction step (if allowed by h)
    if mod(i,h) == 0
        L(:,:,i) = phi(:,:,i+1)*C'*inv(C*phi(:,:,i+1)*C'+V);
        phi(:,:,i+1) = (eye(n_states)-L(:,:,i)*C)*phi(:,:,i+1);
    end
    
    % Update cost sum
    cost_trace(i+1) = cost_trace(i) + trace(P_bar*W)+trace(X_bar*phi(:,:,i+1));
end

% Calculate average cost with trace method
cost = cost_trace(i+1)/(T);
disp(['Average cost = ', num2str(cost)])
end

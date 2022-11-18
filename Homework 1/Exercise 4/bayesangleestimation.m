function [p] = bayesangleestimation(aint,yint,alpha,beta,gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimal control and dynamic programming    %
%Homework 1                                 %
%Assignment 4                               %
%Date: 30-11-2018                           %
%Group: 2                                   %
%Bob Clephas            | 1271431           %
%Tom van de laar        | 1265938           %
%Job Meijer             | 1268155           %
%Marcel van Wensveen    | 1253085           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(alpha);
k = length(yint);

p = zeros(N,k);

for i = 1:k
    for l = 1:N
        theta = l-1;
        g = gamma(mod(theta-yint(i),N)+1);
        if(i == 1)
            a = alpha(l);
            b = 1;
        else
            a = 1;
            b = 0;
            for d = 1:N                
                diff = mod(theta-aint-(d-1),N)+1;
                bd = beta(diff)*p(d,i-1);
                b = b + bd;
            end
        end
        p(l,i) = a*b*g;
    end
    p(:,i) = p(:,i)./sum(p(:,i));
end
end
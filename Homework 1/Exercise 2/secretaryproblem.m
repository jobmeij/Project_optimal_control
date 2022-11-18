function [r] = secretaryproblem(n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimal control and dynamic programming    %
%Homework 1                                 %
%Assignment 2                               %
%Date: 30-11-2018                           %
%Group: 2                                   %
%Bob Clephas            | 1271431           %
%Tom van de laar        | 1265938           %
%Job Meijer             | 1268155           %
%Marcel van Wensveen    | 1253085           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vk = zeros(1,n);

for m = 1:n
    sum = 0;
    for i = m:n-1
        sum = sum + (1/i);
    end
    vk(m) = (m/n)*sum;
end

indexes = find(vk==max(vk));

r = indexes(1)+1;

end
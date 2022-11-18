function [thetaest, X, V, C] = viterbiangleestimation(aint,yint,alpha,beta,gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimal control and dynamic programming    %
%Homework 1                                 %
%Assignment 1                               %
%Date: 30-11-2018                           %
%Group: 2                                   %
%Bob Clephas            | 1271431           %
%Tom van de laar        | 1265938           %
%Job Meijer             | 1268155           %
%Marcel van Wensveen    | 1253085           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(alpha);
h = length(yint);

C = ones(N,h)*NaN;
X = C;
V = C;

for k = h:-1:1
    for i = 1:N
        X(i,k) = (i-1);
        cost = zeros(1,N)*0;
        theta = i-1;
        measurementcost = gamma(mod(yint(h-k+1) - theta,N)+1);
        if(k<h)
            for a_distrubance = find(beta ~= 0)'
                theta_plus = mod(theta - aint - (a_distrubance -1),N);
                prevcost = C(theta_plus+1,k+1);
                disturbancecost = beta(a_distrubance);
                cost(theta_plus+1) = measurementcost * disturbancecost*prevcost;
            end
        end
        
        if(k==h)
            cost(i) = -alpha(i)*measurementcost;
        end
        C(i,k) = min(cost);
        indexes = find((cost == min(cost)) & (cost ~= 0));
        
        z = length(indexes);
        if(z >= size(V,3)+1)
            V(:,:,size(V,3)+1:z) = ones(N,h,z-(size(V,3)))*NaN;
        end
        V(i,k,1:z) = indexes-1;
    end
    C(:,k) = C(:,k)./abs(sum(C(:,k)));
end

startthetas = find(C(:,k) == min(C(:,k)))-1;
for i = 1:length(startthetas)
    thetaestimations(i,h) = startthetas(i);
end
for k = 1:(h-1)
    for j = 1:length(thetaestimations(:,h-k+1))
        V_k = squeeze(V(thetaestimations(j,h-k+1)+1,k,:));
        besttheta = V_k(~isnan(V_k));
        for i = 1:length(besttheta)
            if i == 1
                thetaestimations(j,h-k) = besttheta(i);
            elseif i > 1
                n_vec = length(thetaestimations(:,h-k+1))+1;
                thetaestimations(n_vec,:) = thetaestimations(j,:);
                thetaestimations(n_vec,h-k) = besttheta(i);
            end
        end
    end
end
thetaest = min(thetaestimations,[],1)';
end
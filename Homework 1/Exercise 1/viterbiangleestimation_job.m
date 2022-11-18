function [thetaestimations, X, V, C] = viteriangleestimation(aint, yint, alpha, beta, gamma)
N = length(alpha);
h = length(yint);

C = ones(N,h)*NaN;
X = C;
V = C;

for k = h:-1:1
    for i = 1:N
        X(i,k) = (i-1);
        cost = zeros(1,N);
        if(k<h)
            for a_distrubance = find(beta ~= 0)
                theta = i-1;
                theta_plus = mod(theta - aint - (a_distrubance -1),N);
                theta_diff = mod(theta-theta_plus,N);
                
                prevcost = C(theta_plus+1,k+1);
                diff = mod(yint(h-k+1) - theta,N);
                measurementcost = gamma(diff+1);
                if(measurementcost == 0)
                    %measurementcost = nan;
                end
                disturbancecost = beta(a_distrubance);
                cost(theta_diff+1) = - (measurementcost * disturbancecost) + prevcost;
            end
        end
        
        if(k==h)
            cost(i) = -alpha(i);
            %cost(cost == 0) = nan;
        end
        
        C(i,k) = min(cost);
        indexes = find(cost == min(cost));
        z = length(indexes);
        if(z > size(V,3)+1)
            V(:,:,size(V,3)+1:z) = ones(N,h,z-(size(V,3)))*NaN;
        end
        V(i,k,1:z) = i-X(indexes,k)-1;
    end
end

% bestthetas = find(C(:,k) == min(C(:,k)))-1;
% npaths = length(bestthetas);
% for s = 1:length(bestthetas)
%     for k = 1:h
%         npaths = npaths + length(V(s+1,k));
%     end
% end
% thetaestimations = zeros(length(bestthetas),length(yint));
% for est = 1:length(bestthetas)
%     besttheta = bestthetas(est);
%     for k = 1:h
%         thetaestimations(est,h-k+1) = X(besttheta+1,k);
%         besttheta = mod(besttheta - V(besttheta+1,k,1),N);
%     end
% end
thetaestimations = yint;
end
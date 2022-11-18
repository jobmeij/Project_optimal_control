function draw_diagram(X,C,V)
nstages  = size(X,2);
N = size(X,1);
for k = 1:nstages
    nstates = length(X(~isnan(X(:,k)),k));
    for i = 1:nstates
        if(k < nstages)
            for t = 1:size(V,3)
                quiver(k-1, X(i,k),1,V(i,k,t)-X(i,k),1,'k','MaxHeadsize',0.1);
                %if(~isnan(V(i,k,t)))
                %    annotation('arrow',[k-1 k]./nstages,[X(i,k) V(i,k,t)]./N);
                %end
                hold on;
            end
        end
        plot(k-1,X(i,k),'k.','markersize', 5);
        hold on;
        %text(k-1,X(i,k),num2str(X(i,k)),'HorizontalAlignment','Center');
        text(k-1,X(i,k)+0.4,num2str(C(i,k)),'HorizontalAlignment','Center','color',[1 0 0],'Fontsize',15);
    end
end

xlim([-1, nstages]);
ylim([min(min(X))-1, max(max(X))+1]);
end
function [xix,xiy,xiz] = lqrtrajgeneration(P,T,X0,rho,tau1,N12)
    % Compute constants
    h = floor(T(end)/tau1);
    tau2 = tau1/N12;
    
    r = zeros(3,1:h);
    
    dt = 0.001;
    
    
    
    for i = 1:N-1
        for t = T(i):dt:T(i+1)


        end
    end
    
    
    
    % Generate output
    xix = 
    xiy = 
    xiz = 
end


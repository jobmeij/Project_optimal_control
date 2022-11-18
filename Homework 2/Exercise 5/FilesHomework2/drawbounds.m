function drawbounds(A,B,n,H,Xmin,Xmax,Ymin,Ymax,col)
Ncol = length(col);
for k = 1:H
    for i = 1:n
        Bi = i+(k-1)*n;
        rx = A(i+(k-1)*n,2*k-1);
        ry = A(i+(k-1)*n,2*k);
        if abs(rx) < abs(ry)
            plot([Xmin Xmax],[(B(Bi)-rx*Xmin)/ry (B(Bi)-rx*Xmax)/ry],col(mod(i-1,Ncol)+1));
            hold on;
        else
            plot([(B(Bi)-ry*Ymin)/rx (B(Bi)-ry*Ymax)/rx],[Ymin Ymax],col(mod(i-1,Ncol)+1));
            hold on;
        end
    end
end
end
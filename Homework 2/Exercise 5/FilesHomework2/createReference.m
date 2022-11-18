function ref = createReference(Track, dist)
ref = [];
done = false;
i = 1;
previ = 1;
while done == false
    if isempty(ref)
        point1 = Track(5:6,1);
    else
        point1 = ref(1:2,end);
    end
    
    bestdist = inf;
    bestj = 1;
    tracki = 1;
    for j_ = i:size(Track,2)
        j = mod(j_,size(Track,2))+1;
        testvect = Track(5:6,j)-point1;
        testdist = norm(testvect);
        if testdist < bestdist && j~=i
            bestdist = dist;
            bestj = j;
            tracki = j_;
        end
    end
    i = bestj;
    point2 = Track(5:6,i);
    vect = point2-point1;
    vectlength = norm(vect);
    vect = dist*vect/vectlength;
    
    for d = 1:ceil(vectlength/dist)
        if isempty(ref)
            ref = [ref [point1+vect;tracki]];
        else
            ref = [ref [ref(1:2,end)+vect;tracki]];
        end
    end
    
    if size(ref,2)>=2
        if previ > i
            done = true;
        end
    end
    previ = i;
end
end
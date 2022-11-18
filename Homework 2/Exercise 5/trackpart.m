function part = trackpart(Track, pos, prevpart)
part = prevpart;
for i = 1:size(Track,2)
    j = mod(i,size(Track,2))+1;
    if(inpolygon(pos(1),pos(2),[Track(1,[i j]) Track(3,[j i])],[Track(2,[i j]) Track(4,[j i])]))
        part = i;
        return;
    end
end

end
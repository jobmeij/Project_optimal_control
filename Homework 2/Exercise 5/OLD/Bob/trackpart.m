function part = trackpart(Track, pos, prevpart)
part = prevpart;
for i = 1:size(Track,2)-1
    if(inpolygon(pos(1),pos(2),[Track(1,[i i+1]) Track(3,[i+1 i])],[Track(2,[i i+1]) Track(4,[i+1 i])]))
        part = i;
        return;
    end
end

end
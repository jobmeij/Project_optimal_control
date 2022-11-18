function inside = OnTrack(Track, pos)
    outsideline = inpolygon(pos(1),pos(2),Track(1,:),Track(2,:));
    insideline = ~inpolygon(pos(1),pos(2),Track(3,:),Track(4,:));
    inside = outsideline && insideline;
end
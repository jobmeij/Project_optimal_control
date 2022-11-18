function pass = passfinish(finish,p1,p2)
    pass = false;
    if p1(1) >= finish(1,1) && p1(1) <= finish(1,2) && p2(1) >= finish(1,1) && p2(1) <= finish(1,2)
        if p1(2) >= finish(2,1) && p2(2) <= finish(2,1) 
            pass = true;
            disp("Finish line passed");
        end
    end
end
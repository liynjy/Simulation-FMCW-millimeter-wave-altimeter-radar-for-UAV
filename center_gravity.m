function gx = center_gravity(c)

    gx = 0;
    avg = 0;
    
    N=max(size(c));
    for k=1:N
            gx = gx + k*c(k);
            avg = avg + c(k);
    end
    
    if avg~=0
        gx = gx/avg;
    else
        gx = (N+1)/2;
    end
    
end
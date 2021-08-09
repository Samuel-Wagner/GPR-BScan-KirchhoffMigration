function min_ind = min1(t,M)
    % t is sorted with variation along dim3
    % M is the object to match to
    Np = numel(t);

    lo = 1;
    hi = Np;
    while(lo<hi)
        mid = floor((lo+hi)/2);
        tq = t(mid);
        if(tq==M)
            break;
        elseif(tq > M) %shift down
            hi = mid-1;
        else
            lo = mid+1;
        end
    end

    min_ind = uint32(mid);
end
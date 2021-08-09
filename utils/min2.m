function min_inds = min2(t,M)
    % t is sorted with variation along dim3
    % M is the object to match to
    Nr=numel(M);
    Np = numel(t);
    min_inds = zeros(1,Nr,'uint32');

    for ii = 1:Nr
        Mq = M(ii);
        lo = 1;
        hi = Np;
        while(lo<hi)
            mid = ceil((lo+hi)/2);
            tq = t(mid);
            if(tq==Mq)
                break;
            elseif(tq > Mq) %shift down
                hi = mid-1;
            else
                lo = mid+1;
            end
        end
        min_inds(ii) = uint32(mid);
    end
end
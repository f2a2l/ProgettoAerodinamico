function [mind] = getMinDist(x1,y1,x2,y2)

    %% user settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minX = 0.8;
    maxX = 1.2;
    
    
    %% reduce vectors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx1 = find(x1 > minX);
    idx2 = find(x2 <= maxX);
    x1 = x1(idx1);
    y1 = y1(idx1);
    x2 = x2(idx2);
    y2 = y2(idx2);
    
    %% find min distance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = length(x1);
    n = length(x2);
    mind = norm([(x1(1)-x2(1)), (y1(1)-y2(1))]);
    for ii = 1:m
        for jj = 1:n
            local_d = norm([(x1(ii)-x2(jj)), (y1(ii)-y2(jj))]);
            if local_d < mind
                mind = local_d;
            end
        end
    end
    
end


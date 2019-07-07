function [y_out] = interpFun(t,y,required_t)
%interpFun Returns value of function at time required_t by linear
%interpolation; t must be a vector, as well as y

    if ~isvector(t) || ~isvector(y) || length(t) ~= length(y)
        error('last two inputs must be vectors of equal size.')
    end
    
    m = length(t);
    
    ii=find(t >= required_t);

    if isempty(ii)

        y_out = y(end);
        if required_t ~= t(end)
            warning('out of bounds; returning closest value.')
        end
        
    else

        ii=ii(1);

        if ii == 1 
            y_out = y(1);
            if required_t ~= t(1)
                warning('out of bounds; returning closest value.')
            end
        else
            m = (required_t - t(ii-1)) / (t(ii) - t(ii-1));
            y_out = y(ii-1) + m*(y(ii)-y(ii-1));
        end

    end

end
    
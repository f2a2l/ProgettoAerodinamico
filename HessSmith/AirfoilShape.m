function [xR, yR] = AirfoilShape (aname, varargin)
%


if (~isempty (varargin))
    method = varargin{1};
    if (length(varargin)>1)
        n = varargin{2};
    else
        n = 49;
    end
else
    method = 'constant';
    n      = 49;
end

% Naca identification
type = aname(5:end);
        
% 4 digits
if (length(type) < 5)
    
    maxord = str2num (type(1));
    posmax = str2num (type(2));
    thich  = str2num (type(3:end));
    [x1, y1, x2, y2] = NACA_4d (maxord, posmax, thich, n, method);
    
% 5 digits
else
    
    series = str2num (type(1:3));
    thich  = str2num (type(4:end));
    [x1, y1, x2, y2] = NACA_5d (series, thich, n, method);

end

xd = x1(:,1);
yd = y1(:,1);
xv = x2(:,1);
yv = y2(:,1);

xR = flipud([flipud(xd); xv(2:end)]);
yR = flipud([flipud(yd); yv(2:end)]);

return
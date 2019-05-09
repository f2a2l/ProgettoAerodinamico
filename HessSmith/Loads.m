function [Cl, Cd, Cm, varargout] = Loads (p, Cp, alpha, varargin)


n = length(p.panel);
if (~isempty(varargin))
dist = varargin{1};
else
    dist = [0 0];
end
CFx = 0;
CFy = 0;
Cm = 0;
CmLE = 0;
xref = dist(1) + 0.25;
yref = 0;
for i = 1:n-1
    
    dx = p.panel(i).P2.x - p.panel(i).P1.x;
    dy = p.panel(i).P2.y - p.panel(i).P1.y;
    
    CFx = CFx + Cp(i)*dy;
    CFy = CFy - Cp(i)*dx;
    %Cm  = Cm + Cp(i)*(dx*(p.panel(i).C.x + dy*p.panel(i).C.y));
    Cm  = Cm + Cp(i)*(dx*((p.panel(i).C.x - xref) + dy*(p.panel(i).C.y - yref)));
    CmLE  = CmLE + Cp(i)*(dx*(p.panel(i).C.x - dist(1)) + dy*(p.panel(i).C.y- dist(2)));
    
end
varargout{1} = CmLE;
Cd = [CFx CFy]*[cos(deg2rad(alpha)); sin(deg2rad(alpha))];
Cl = [CFx CFy]*[-sin(deg2rad(alpha)); cos(deg2rad(alpha))];

return
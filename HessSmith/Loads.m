function [Cl, Cd] = Loads(p, Cp, alpha)


n = length(p.panel);

CFx = 0;
CFy = 0;

for i = 1:n-1
    
    dx = p.panel(i).P2.x - p.panel(i).P1.x;
    dy = p.panel(i).P2.y - p.panel(i).P1.y;
    
    CFx = CFx + Cp(i)*dy;
    CFy = CFy - Cp(i)*dx;
    
end
Cd = [CFx CFy]*[cos(deg2rad(alpha)); sin(deg2rad(alpha))];
Cl = [CFx CFy]*[-sin(deg2rad(alpha)); cos(deg2rad(alpha))];

return
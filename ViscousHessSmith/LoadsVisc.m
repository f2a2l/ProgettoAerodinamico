function Cd = LoadsVisc(p, Cf, alpha)

    
n = length(p.panel);

CFx = 0;
CFy = 0;

for i = 1:n-1
    
    dx = abs(p.panel(i).P2.x - p.panel(i).P1.x);    
    CFx = CFx + Cf(i)*dx;

end

Cd = [CFx CFy]*[cos(deg2rad(alpha)); sin(deg2rad(alpha))];

return
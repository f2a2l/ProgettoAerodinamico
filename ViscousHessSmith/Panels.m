function [p] = Panels(x, y)
%Build the data structure containing all the data related to the panels:
%end points location, center point location, panel length, angle and 
%rotation matrix. Adds also a wake panel.


n = length(x);

p = struct;

for i = 1:n-1
  
    p.panel(i).P1.x  = x(i);
    p.panel(i).P2.x  = x(i+1);
    p.panel(i).P1.y  = y(i);
    p.panel(i).P2.y  = y(i+1);
    p.panel(i).C.x   = (x(i)+x(i+1))/2;
    p.panel(i).C.y   = (y(i)+y(i+1))/2;
    p.panel(i).d     = sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2);
    p.panel(i).beta  = atan2 ((y(i+1)-y(i)), (x(i+1)-x(i))); 
    p.panel(i).R     = Rotation(p.panel(i).beta);
    
end

[p] = addWakePanel (p);

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [p] = addWakePanel(p)

n = length(p.panel);

panel_1 = p.panel(1);
panel_2 = p.panel(n);

p.panel(n+1).P1.x = panel_1.P1.x;
p.panel(n+1).P1.y = panel_1.P1.y;
p.panel(n+1).P2.x = (panel_1.P1.x + (panel_1.P1.x - panel_1.P2.x) +     ...
                    panel_2.P2.x + (panel_2.P2.x - panel_2.P1.x))/2;
p.panel(n+1).P2.y = (panel_1.P1.y + (panel_1.P1.y - panel_1.P2.y) +     ...
                    panel_2.P2.y + (panel_2.P2.y - panel_2.P1.y))/2;

p.panel(n+1).d    = sqrt((p.panel(n+1).P2.x-p.panel(n+1).P1.x)^2 +      ...
                        (p.panel(n+1).P2.y-p.panel(n+1).P1.y)^2);
p.panel(n+1).C.x  = (p.panel(n+1).P1.x+p.panel(n+1).P2.x)/2;
p.panel(n+1).C.y  = (p.panel(n+1).P1.y+p.panel(n+1).P2.y)/2;
p.panel(n+1).beta = atan2 ((p.panel(n+1).P2.y-p.panel(n+1).P1.y),       ...
                           (p.panel(n+1).P2.x-p.panel(n+1).P1.x));
p.panel(n+1).R     = Rotation(p.panel(n+1).beta);                         

return
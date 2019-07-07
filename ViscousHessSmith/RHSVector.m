function [RHS] = RHSVector (p, alpha, U)


n = length(p.panel);

uInf = U*[cos(deg2rad(alpha)); sin(deg2rad(alpha))];

RHS = zeros (n, 1);

for i = 1:n
    RHS(i,1) = -[-sin(p.panel(i).beta) cos(p.panel(i).beta)]*uInf;
end

return
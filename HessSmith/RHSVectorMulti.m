function [RHS] = RHSVectorMulti (p, alpha, U)


n = length(p.panel);

alpha1 = alpha(1); % il sistema di riferimento globale è solidale al primo profilo!

uInf1 = U*[cos(deg2rad(alpha1)); sin(deg2rad(alpha1))];


RHS = zeros (n, 1);

for i = 1:n
   
    RHS(i,1) = -[-sin(p.panel(i).beta) cos(p.panel(i).beta)]*uInf1;
   
end

return
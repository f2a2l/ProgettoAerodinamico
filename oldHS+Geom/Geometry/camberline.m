
function [xc,yc] = camberline(c1,c2,c3,c4,k)

    xc = (3*c1*k*(1-k)^2)+(3*c2*(1-k)*k^2)+k^3;
    yc = (3*c3*k*(1-k)^2)+(3*c4*(1-k)*k^2);

end 
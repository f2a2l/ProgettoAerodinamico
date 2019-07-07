function [AIC] = AICMatrix (p)




n = length(p.panel);

AIC = zeros (n, n);

for i = 1:n
    
    ni = [-sin(p.panel(i).beta) cos(p.panel(i).beta)];
    
    for j = 1:n-1
        
        if (i == j)
            
            Rt = Rotation(p.panel(j).beta)';
            us = [0; 1/2];
            uv = [1/2; 0];
            us = Rt*us;
            uv = Rt*uv;
            
        else
            
            us = ConstantSource2D (1, p.panel(j), p.panel(i));
            uv = [us(2); -us(1)];
            
        end
        
        AIC(i,j) = ni*us;
        AIC(i,n) = AIC(i,n) + ni*uv;
        
    end
end

return
function [AIC] = AICMatrixMulti (p, nairfoils)


ntot = length(p.panel);

AIC = zeros(ntot, ntot);

nterz = (ntot - nairfoils)/nairfoils;

for i = 1:ntot


    ni = [-sin(p.panel(i).beta) cos(p.panel(i).beta)];


    % sotto matrice Aij (sorgenti su sorgenti) e bs (sorgenti su vortici)
    for j = 1:ntot-nairfoils

        if (i == j)

            Rt = Rotation(p.panel(j).beta)';
            us = [0; 1/2];
            us = Rt*us;

        else
    
            us = ConstantSource2D(1, p.panel(j), p.panel(i));

        end

        AIC(i,j) = ni*us;

    end


    %sottomatrice b_v (vortici su sorgenti) 
    for k = 1:nairfoils

        for j = ((k-1)*nterz + 1 ):k*nterz
            if (i == j)

                Rt = Rotation(p.panel(j).beta)';
                uv = [1/2; 0];
                uv = Rt*uv;

            else

                us = ConstantSource2D(1, p.panel(j), p.panel(i));
                uv = [us(2); -us(1)];

            end

            AIC(i, ntot-nairfoils + k) = AIC(i, ntot- nairfoils + k) + ni*uv;

        end
    end


end


return
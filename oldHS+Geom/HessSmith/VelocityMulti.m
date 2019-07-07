function[vaux] = VelocityMulti (p, p1, alpha, U, SOL)


ntot = length(p.panel);
nairfoils = length(p1);
nterz = (ntot - nairfoils)/nairfoils;

alpha1 = alpha(1); % il sist. di riferimento � solidale al primo profilo
uInf1 = U*[cos(deg2rad(alpha1)); sin(deg2rad(alpha1))];

v = zeros(ntot-nairfoils,1);

for i = 1:ntot-nairfoils 

    ti = [cos(p.panel(i).beta) sin(p.panel(i).beta)];

    %CONTRIBUTO SORGENTI 
    
    for j = 1:ntot-nairfoils


        if (i == j)

            Rt = Rotation(p.panel(j).beta)';
            us = Rt*[0; 1/2]*SOL(j);
            % uv = Rt*[1/2; 0]*SOL(end);

        else

            us = ConstantSource2D (SOL(j), p.panel(j), p.panel(i));


        end

        v(i) = v(i) + ti*us ; %+ ti*uv;

    end
    
    % CONTRIBUTO VORTICI
    
    for k = 1:nairfoils


        for j =  ((k-1)*nterz + 1 ):k*nterz

            if (i == j)

                Rt = Rotation(p.panel(j).beta)';
                uv = Rt*[1/2; 0]*SOL(end-nairfoils+k);

            else

                uv = ConstantVortex2D (SOL(end-nairfoils+k), p.panel(j), p.panel(i));


            end

            v(i) = v(i) + ti*uv;


        end
    end


   
    v(i) = v(i) + ti*uInf1;
   

end


vaux = zeros(nterz,nairfoils);

for k = 1:nairfoils
    j = 0;
    for  i = ((k-1)*nterz + 1):k*nterz
        j = j + 1;
        vaux(j,k) = v(i);
    end
end





return


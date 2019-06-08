function [vaux] = VelocityMulti (p, metaPan, nairfoils, alpha, U, SOL)

% extract input
npan = metaPan.npan;
idx_zeroPan = metaPan.idx_zeroPan;

ntot = length(p.panel);

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


        for j =  (1:npan(k)) + idx_zeroPan(k)

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


vaux = cell(1,nairfoils);

for k = 1:nairfoils
    vaux{k} = v((1:npan(k)) + idx_zeroPan(k));
end



return


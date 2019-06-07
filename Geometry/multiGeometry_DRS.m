function [x, y, totLength] = multiGeometry_DRS(npoint, aname, alpha, dist, crel, varargin)

    nairfoils = length(alpha);
    
    dist(:,1) = dist(:,1) + 1;

    if isempty(varargin)
        pltFlag = false;
    else
        pltFlag = varargin{1};
    end

    % preallocation of airfoil shape
    x = zeros(2*npoint - 1, nairfoils);
    y = zeros(2*npoint - 1, nairfoils);

    % get airfoil shape
    for i = 1:nairfoils
        [x(:,i), y(:,i)] = AirfoilShape(aname(i,:), npoint);
        if i > 1
            x(:,i) = crel(i-1) * x(:,i);
            y(:,i) = crel(i-1) * y(:,i);
        end
    end

    % resize airfoils and apply angle of attack
    for i = 2:nairfoils

        alpha_rel = alpha(i) - alpha(1);

        R= [cos(alpha_rel*pi/180) sin(alpha_rel*pi/180);
            -sin(alpha_rel*pi/180) cos(alpha_rel*pi/180)];

        coord_mat = [x(:,i), y(:,i)]';
        coord_mat = R * coord_mat;
        coord_mat = coord_mat';

        x(:,i) = coord_mat(:,1) + dist(i-1,1);
        y(:,i) = coord_mat(:,2) + dist(i-1,2);

    end

    
    DRS_angle = alpha(2) - alpha(1);
    [x, y] = DRS_OPEN(x,y,DRS_angle);
    
    % plot airfoils
    if pltFlag
        figure('Name','Multi element profile geometry','NumberTitle','off')
        hold on
        axis equal
        for i = 1:nairfoils
            scatter(x(:,i),y(:,i), '.', 'black')
        end
        hold off
    end

    % calculate total length (normalised on first airfoil's parameter)
    totLength = max(x(:));

end
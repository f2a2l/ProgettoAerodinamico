function [x, y, totLength, yTotLength] = multiGeometry(npoint, aname, alpha, dist, crel, varargin)

    nairfoils = length(alpha);

    dist(:,1) = dist(:,1) + 1;

    if isempty(varargin)
        pltFlag = false;
    else
        pltFlag = varargin{1};
    end

    % preallocation of airfoil shape
    x = cell(1,2);
    y = cell(1,2);

    orig_npoint = npoint;

    % get airfoil shape
    for i = 1:nairfoils

        % adjust number of points
        if i > 1
            npoint = getNpt(orig_npoint, crel(i-1));
        end

        % get shape
        [xcurr, ycurr] = AirfoilShape(aname(i,:), npoint);

        % apply geometric transformation to rear profiles
        if i > 1

            % adjust airfoil size
            xcurr = crel(i-1) * xcurr;
            ycurr = crel(i-1) * ycurr;

            % rotation
            alpha_rel = alpha(i) - alpha(1);
            R= [cos(alpha_rel*pi/180) sin(alpha_rel*pi/180);
                -sin(alpha_rel*pi/180) cos(alpha_rel*pi/180)];
            coord_mat = [xcurr, ycurr]';
            coord_mat = R * coord_mat;
            coord_mat = coord_mat';

            % translation
            xcurr = coord_mat(:,1) + dist(i-1,1);
            ycurr = coord_mat(:,2) + dist(i-1,2);

        end

        % save
        x{i} = xcurr;
        y{i} = ycurr;

    end
    
    % plot airfoils
    if pltFlag
        figure('Name','Multi element profile geometry','NumberTitle','off')
        hold on
        axis equal
        for i = 1:nairfoils
            scatter(x{i},y{i}, '.', 'black')
        end
        hold off
    end

    % calculate total length (normalised on first airfoil's parameter)
    totLength = 0;
    for ii = 1:length(x)
        totLength = max(totLength, max(x{ii}));
    end
    for ii = 1:length(x)
        a  = x{ii};
        idx = find(a == totLength, 1);
        if ~isempty(idx)
            yTotLength = a(idx);
        end
    end

end

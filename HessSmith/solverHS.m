function [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, varargin)
% Usage:
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha)
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, pltFlag)
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, dist, crel)
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, dist, crel, pltFlag)

% maxdCp is a matrix; each row corresponds to an airfoil, first column corresponds to
% lower part, second column corresponds to upper part

% TODO: input check => mi conviene davvero farlo? 
% magari lo faccio e poi lo commento prima di inserirlo in ciclo di ottimizzatore


if length(alpha) == 1

    if isempty(varargin)
        pltFlag = false; % default: do not plot
    else
        pltFlag = varargin{1};
    end

    aname1 = aname(1, :);
    alpha1 = alpha;

    % Airfoil discretization and plotting
    [x, y] = AirfoilShape(aname1, npoint);
    [p1] = Panels(x, y);
    if pltFlag
        figure('Name','Profile geometry','NumberTitle','off')
        plot(x,y, 'black')
        axis equal
        hold off
    end

    % Aerodynamic Influence Coefficients Matrix [AIC]
    [AIC] = AICMatrix (p1);

    % Right Hand Side Vector {RHS}
    [RHS] = RHSVector(p1, alpha1, 1);

    % System solution
    SOL = AIC\RHS;

    % Velocity
    [v] = Velocity(p1, alpha1, 1, SOL); % need this for Cp

    % Pressure - FIXME: fix calculation of maxdCp after reading Valarezo-Chin
    maxdCp = zeros(1,2);
    [Cp] = PressureCoeff(v, 1);
    ncp = floor(size(Cp,1)/2);
    maxdCp(1,1) = max(Cp(1:ncp,1)) - min(Cp(1:ncp,1));
    ncp = ncp + 1;
    maxdCp(1,2) = max(Cp(ncp:end,1)) - min(Cp(ncp:end,1));



    % Aerodynamic coefficients
    [Cl, Cd] = Loads(p1, Cp, alpha1); % omitted arguments: Cm, CmLE

    
elseif (isempty(varargin)) && length(alpha) > 1
    error('no distance specified!')


elseif (~isempty(varargin)) && length(alpha) >= 2

    % Identify case
    nairfoils = length(alpha);

    % Load inputs
    if length(varargin) == 2
        dist = varargin{1};
        crel = varargin{2};
        pltFlag = false; % default: do not plot
    elseif length(varargin) == 3
        dist = varargin{1};
        crel = varargin{2};
        pltFlag = varargin{3}; % default: do not plot
    end
    
    % security check FIXME: sicuramente sbagliato, forse togli
    for i = 1 : nairfoils-2
        if dist(i+1,1) - dist(i,1) < 1 && abs(dist(i,2) - dist(i+1,2)) < 0.5
            error('Airfoils are too close, overlapped or not in the right order!')
        end
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

    % plot airfoils
    if pltFlag
        figure('Name','Multi element profile geometry','NumberTitle','off')
        hold on
        axis equal
        for i = 1:nairfoils
            plot(x(:,i),y(:,i), 'black')
        end
        hold off
    end

    % panels
    for i = nairfoils:-1:1 % this runs backwards to avoid preallocation issues!
        [p1(i)] = Panels(x(:,i), y(:,i));
    end
    [p] = PanelsMulti(p1);

    % Influence matrix
    [AIC] = AICMatrixMulti(p, nairfoils);

    % Right Hand Side Vector {RHS}
    [RHS] = RHSVectorMulti (p, alpha, 1);

    % Solution
    SOL = AIC\RHS;

    % Velocity on profile
    [v] = VelocityMulti(p, p1, alpha, 1, SOL);

    % Preallocation of coefficients
    [np, nairfoils] = size(v);
    Cp = zeros(np, nairfoils);
    maxdCp = zeros(nairfoils, 2);
    Cl = zeros(nairfoils,1);
    Cd = zeros(nairfoils,1);

    % Calculation of coefficients
    for i = 1:nairfoils
        [Cp(:,i)] = PressureCoeff(v(:,i), 1);
        if i == 1
            [Cl(i), Cd(i)] = Loads(p1(i), Cp(:,i), alpha(i));
        else
            [Cl(i), Cd(i)] = Loads(p1(i), Cp(:,i), alpha(i));
        end
        % FIXME: fix calculation of maxdCp after reading Valarezo-Chin
        ncp = floor(size(Cp,1)/2);
        maxdCp(i, 1) = max(Cp(1:ncp, i)) - min(Cp(1:ncp, i));
        ncp = ncp + 1;
        maxdCp(i, 2) = max(Cp(ncp:end, i)) - min(Cp(ncp:end, i));
    end


else  
    error('ERROR! Wrong input number')
end


return

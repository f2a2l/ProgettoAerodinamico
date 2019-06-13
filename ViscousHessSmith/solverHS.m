function [Cl, Cd, xmax, ymax, Cp, v, maxdCp, x, y, p, p1, SOL, metaPan, nairfoils] = solverHS(npoint, aname, alpha, varargin)
% Usage:
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha)
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, dist, crel)

% maxdCp is a matrix; each row corresponds to an airfoil, first column corresponds to
% lower part, second column corresponds to upper part

% TODO: input check => mi conviene davvero farlo? 
% magari lo faccio e poi lo commento prima di inserirlo in ciclo di ottimizzatore

metaPan = 0;

if length(alpha) == 1 && isempty(varargin)

    p1 = 0;

    nairfoils = 1;
    xmax = 1;
    ymax = 0;

    aname1 = aname(1, :);
    alpha1 = alpha;

    % Airfoil discretization and plotting
    [x, y] = AirfoilShape(aname1, npoint);
    [p] = Panels(x, y);

    % Aerodynamic Influence Coefficients Matrix [AIC]
    [AIC] = AICMatrix (p);

    % Right Hand Side Vector {RHS}
    [RHS] = RHSVector(p, alpha1, 1);

    % System solution
    SOL = AIC\RHS;

    % Velocity
    [v] = Velocity(p, alpha1, 1, SOL); % need this for Cp

    % Pressure
    [Cp] = PressureCoeff(v, 1);
    Cpmin = min(Cp);
    maxdCp = abs(Cpmin - Cp(end));



    % Aerodynamic coefficients
    [Cl, Cd] = Loads(p, Cp, alpha1); % omitted arguments: Cm, CmLE

    % Stuff to achieve uniformity with multiple airfoil case
    x = {x}; y = {y};
    Cp = {Cp}; v = {v};


elseif (~isempty(varargin)) && length(alpha) >= 2

    % Identify case
    nairfoils = length(alpha);

    % Load inputs
    dist = varargin{1};
    crel = varargin{2};

    % get multi geometry
    [x,y, xmax, ymax] = multiGeometry(npoint, aname, alpha, dist, crel);

    % panels
    for i = 1:nairfoils % this runs backwards to avoid preallocation issues!
        [p1(i)] = Panels(x{i}, y{i});
    end
    [p, metaPan] = PanelsMulti(p1);

    % Influence matrix
    [AIC] = AICMatrixMulti(p, metaPan, nairfoils);

    % Right Hand Side Vector {RHS}
    [RHS] = RHSVectorMulti(p, alpha, 1);

    % Solution
    SOL = AIC\RHS;

    % Velocity on profile
    v = VelocityMulti(p, metaPan, nairfoils, alpha, 1, SOL);

    % Preallocation of coefficients
    Cp = cell(1, nairfoils);
    maxdCp = zeros(nairfoils, 1);
    Cl = zeros(nairfoils,1);
    Cd = zeros(nairfoils,1);

    % Calculation of coefficients
    for i = 1:nairfoils
        Cp{i} = PressureCoeff(v{i}, 1);
        [Cl(i), Cd(i)] = Loads(p1(i), Cp{i}, alpha(i));
        % calculation of maxdCp for Valarezo-Chin
        c = Cp{i};
        minCp = min(c);
        maxdCp(i) = abs(minCp-c(end));
    end
    

else  
    error('wrong input; see documentation for instructions on how to use this function.')
end


return

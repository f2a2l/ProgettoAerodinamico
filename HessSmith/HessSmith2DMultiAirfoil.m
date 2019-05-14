function [x, y, p1, Cp , v, varargout] = HessSmith2DMultiAirfoil(npoint, U, aname, alpha, varargin)% dist, crel)

t1 = cputime;


if (isempty (varargin)) && length(alpha) == 1

    aname1 = aname(1, :);
    fprintf('Single airfoil case')
    t1 = cputime;

    alpha1 = alpha;
    % Airfoil discretization
    [x, y] = AirfoilShape(aname1, npoint);
    [p1] = Panels (x, y);

    % Aerodynamic Influence Coefficients Matrix [AIC]
    [AIC] = AICMatrix (p1);

    % Right Hand Side Vector {RHS}
    [RHS] = RHSVector (p1, alpha1, U);

    % System solution
    SOL = AIC\RHS;

    % Velocity
    [v] = Velocity (p1, alpha1, U, SOL);

    % Pressure
    [Cp] = PressureCoeff (v, U);

    % Aerodynamic coefficients
    [Cl, Cd, Cm, CmLE] = Loads (p1, Cp, alpha1);
    %CmLE = Cm;
    t2 = cputime;

    % Print
    PrintFun (int2str(1), t2-t1, Cl, Cd, Cm, CmLE);
    
    
     t3 = cputime;
    [uField,vField, Xgrid, Ygrid] = velocityField(p1, p1, alpha, U, SOL, x, y);
    t4 = cputime;
    varargout{2} = uField;
    varargout{3} = vField;
    varargout{4} = Xgrid;
    varargout{5} = Ygrid;
    fprintf('Time for velocity field computation: %6.2f s \n', t4 - t3)

elseif (isempty (varargin)) && length(alpha) > 1
    error('no distance specified!')


elseif (~isempty (varargin)) && length(alpha) >= 2
    % identify case
    nairfoils = length(alpha);
    fprintf('Multi airfoil case: n airfoils is: %6d \n', nairfoils)

    dist = varargin{1};
    c = 1;
    % security check FIXME: sicuramente sbagliato, forse togli
    for i = 1 : nairfoils-2
        if dist(i+1,1) - dist(i,1) < 1 && abs(dist(i,2) - dist(i+1,2)) < 0.5
            error('Airfoils are too close, overlapped or not in the right order!')
        end
    end



    crel = varargin{2};



    x = zeros(2*npoint - 1,nairfoils);
    y = zeros(2*npoint - 1,nairfoils);


    for i = 1:nairfoils
        [x(:,i), y(:,i)] = AirfoilShape (aname(i,:), npoint);
        if i > 1
            x(:,i) = crel(i-1) * x(:,i);
            y(:,i) = crel(i-1) * y(:,i);
        end
    end

    for i = 2:nairfoils
    % TODO: qui ho fatto modifiche, funzionerà davvero?

        alpha_rel = alpha(i) - alpha(1);

        R= [cos(alpha_rel*pi/180) sin(alpha_rel*pi/180);
            -sin(alpha_rel*pi/180) cos(alpha_rel*pi/180)];

        coord_mat = [x(:,i), y(:,i)]';
        coord_mat = R * coord_mat;
        coord_mat = coord_mat';

        x(:,i) = coord_mat(:,1) + dist(i-1,1);
        y(:,i) = coord_mat(:,2) + dist(i-1,2);

    end


    for i = 1:nairfoils
        [p1(i)] = Panels(x(:,i), y(:,i));
    end

    [p] = PanelsMulti(p1);

    varargout{1} = p;
    %Influence matrix
    [AIC] = AICMatrixMulti (p, nairfoils);


    % Right Hand Side Vector {RHS}
    [RHS] = RHSVectorMulti (p, alpha, U);
    SOL = AIC\RHS;
    %Velocity
    [v] = VelocityMulti (p, p1, alpha, U, SOL);

    [np, nairfoils] = size(v);
    Cp = zeros(np, nairfoils);
    Cl = zeros(nairfoils,1);
    Cd = zeros(nairfoils,1);
    Cm = zeros(nairfoils,1);
    CmLE = zeros(nairfoils,1);


    for i = 1:nairfoils
        [Cp(:,i)] = PressureCoeff (v(:,i), U);
        if i == 1
            [Cl(i), Cd(i), Cm(i)] = Loads (p1(i), Cp(:,i), alpha(i));
            CmLE(i) = Cm(i);
        else
            [Cl(i), Cd(i), Cm(i), CmLE(i)] = Loads (p1(i), Cp(:,i), alpha(i), dist(i-1,:));
        end
        t2 = cputime;
        PrintFun (int2str(i), t2-t1, Cl(i), Cd(i), Cm(i), CmLE(i));
    end



    t3 = cputime;
    [uField,vField, Xgrid, Ygrid] = velocityField(p, p1, alpha, U, SOL, x, y);
    t4 = cputime;
    varargout{2} = uField;
    varargout{3} = vField;
    varargout{4} = Xgrid;
    varargout{5} = Ygrid;
    fprintf('Time for velocity field computation: %6.2f s \n', t4 - t3)
    [cPField] = pressureField(uField, vField, U);
    varargout{6} = cPField;
else  error('ERROR! Wrong input number')

end

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [] = PrintFun (aname, time, Cl, Cd, Cm, CmLE)

fprintf ('\n\n  ==== Hess-Smith method ====');
fprintf ('\n  ====     2D Steady     ==== \n\n');
fprintf ('  Airfoil %s \n', aname);
fprintf ('\n  Exec time: %f', time);
fprintf ('\n  Cl = %f', Cl);
fprintf ('\n  Cd = %f', Cd);
fprintf ('\n  Cm = %f', Cm);
fprintf ('\n  CmLE = %f \n\n', CmLE);

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [] = PlotFun (x, y, p, Cp)

figure
hold on
grid on
plot (x, y, '-*r', 'linewidth', 2)
for i = 1:length(p.panel)
    plot (p.panel(i).C.x, p.panel(i).C.y, 'bd', 'linewidth', 2)
end
axis 'equal'

figure
hold on
grid on
plot (x, y, '-r', 'linewidth', 2)
for i = 1:length(p.panel)
    if (i <= ceil(length(p.panel)/2))
        plot (p.panel(i).C.x, -Cp(i), 'bd', 'linewidth', 2)
    else
        plot (p.panel(i).C.x, -Cp(i), 'ko', 'linewidth', 2)
    end
end
axis ([0 1 -1.5 1])

return
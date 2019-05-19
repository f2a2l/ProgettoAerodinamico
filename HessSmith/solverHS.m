<<<<<<< HEAD
function [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, varargin)
% Usage:
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha)
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, pltFlag)
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, dist, crel)
% - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, dist, crel, pltFlag)

% maxdCp is a matrix; each row corresponds to an airfoil, first column corresponds to
% lower part, second column corresponds to upper part

% For testing, use npoint = 49

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
        [x(:,i), y(:,i)] = AirfoilShape (aname(i,:), npoint);
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
=======
% Hess-Smith panel method for multi-airfoil case: the software allows to consider
% single or multi airfoil case, with an indefinite number of airfoils. For
% each airfoil, the name, the angle of attack and the leading edge
% coordinates (in a reference system defined by the body axes of the first airfoil
% and with the origin at its leading edge) must be specified (interactive
% input).
% Please remember: when inserting airfoil name, leave a blank space
% at the end for 4 digits NACA airfoil.
% author: < Federico Messanelli: federico.messanelli@polimi.it >



function solverHS(AirfNumb, aname, alpha, dist, crel)

    if nargin == 3
        dist = 0;
    elseif nargin == 4
        crel = 1;
    end

    % TODO: check inputs size

    npoint = 80;
    U = 1;
    % -------------------------------------------------------------------------
    % Hess-Smith method (external function)

    if AirfNumb == 1
        [x, y, p1, Cp , v, p, uField, vField, Xgrid, Ygrid] = HessSmith2DMultiAirfoil (npoint, U, aname, alpha);
        dist = [0 0];
    else
        [x, y, p1, Cp , v, p, uField, vField, Xgrid, Ygrid] = HessSmith2DMultiAirfoil (npoint, U, aname, alpha, dist, crel);

    end
    

    
    nairfoils = AirfNumb;
    yVel = (min(dist(:,2))-1:0.25:max(dist(:,2))+1)';
    xVel = -1.1*ones(length(yVel),1);
    unitVect = ones(length(yVel),1);
    for j = 1:nairfoils

    %     figure(j)
    %     plot(x(:,j),y(:,j),'b-o')
    %     axis equal
    %     hold on
    %     for i = 1:length(p1(j).panel)
    %         RR=[cos(p1(j).panel(i).beta) sin(p1(j).panel(i).beta); - sin(p1(j).panel(i).beta) cos(p1(j).panel(i).beta) ];
    %         xc =[p1(j).panel(i).C.x
    %             p1(j).panel(i).C.x + RR(2,1)*0.5];
    % 
    %         yc = [p1(j).panel(i).C.y
    %             p1(j).panel(i).C.y + RR(2,2)*0.5];
    %         plot(xc , yc,'k')
    %         hold on
    %     end
        figure(99)
        plot(x(:,j),y(:,j),'-','color',rand(1,3),'linewidth',2)
        axis equal
        axis([-1 max(dist(:,1))+1 min(dist(:,2))-1 max(dist(:,2))+1])
        hold on
        xlabel('x/c')
        ylabel('y/c')
        quiver(xVel, yVel, cos(alpha(1)*pi/180)*unitVect, sin(alpha(1)*pi/180)*unitVect,'r','linewidth',2)

    end
    fileID = fopen('Cp.txt','w');
    for k = 1:nairfoils
    numUp1 = [];
    numLw1 = [];
    n = length(p1(k).panel)-1;

    fprintf(fileID,'Airfoil %9s\r\n',aname(k,:));
    fprintf(fileID,'%3s %3s\r\n','x','Cp');


    for i = 1:n
        if (i <= ceil(n/2))
            numLw1 = [numLw1; p1(k).panel(i).C.x Cp(i,k)];
        else
            numUp1 = [numUp1; p1(k).panel(i).C.x Cp(i,k)];
        end
    end
    numTot = [numLw1;numUp1];
    n2 = length(numTot);
    for i = 1:n2
    fprintf(fileID,'%3.4f %3.4f\r\n',numTot(i,:));
    end
    if min(numUp1(:,2)) < min(numLw1(:,2))
        par = min(numUp1(:,2));
    else
        par =  min(numLw1(:,2));
    end
    figure(100+k)
    hold on
    grid on
    plot (numUp1(:,1), numUp1(:,2), '-r', 'linewidth', 2)
    hold on
    plot (numLw1(:,1), numLw1(:,2), '-b', 'linewidth', 2)
    %plot (x, y, '-k', 'linewidth', 2)
    % legend ( 'Num Upper', 'Num Lower')
    str = sprintf('Airfoil %d, alpha = %d',k,alpha(k));
    title (str)
    xlabel ('c')
    ylabel ('Cp')
    if k == 1
        axis([0  1 min(numUp1(:,2)) 1])
    else
    axis([0 + dist(k-1,1) 1+dist(k-1,1) par 1])
    end
    legend ( 'Upper Side', 'Lower Side')
    fprintf(fileID,'\r\n');
    end
    fclose(fileID);
    % if alpha(1) == 0
    starty = [-1.5:0.025:1.5]';
    startx = -0.48*ones(length(starty),1);
    % elseif alpha(1) < 0
    %     
    %     startx = [-1:0.05:1]';
    %     starty = 1*ones(length(startx),1);
    % else
    %      startx = [-1:0.05:1]';
    %     starty = -1*ones(length(startx),1);
    % end
    %XY = stream2(Xgrid,Ygrid,uField,vField,startx,starty);
    % figure()
    % hlines=streamline(Xgrid,Ygrid,uField,vField,startx,starty);
    % set(hlines,'LineWidth',2,'Color','r')
    % hold on
    % 
    % plot(Xgrid, Ygrid,'k.')
    % 
    % for j = 1:nairfoils
    % 
    %     plot(x(:,j),y(:,j),'-','color',rand(1,3),'linewidth',2)
    %     axis equal
    %     axis([-1 dist(end,1)+1.5 min(dist(:,2))-0.5 max(dist(:,2))+0.5])
    %     hold on
    %     xlabel('x/c')
    %     ylabel('y/c')
    %     quiver(xVel, yVel, cos(alpha(1)*pi/180)*unitVect, sin(alpha(1)*pi/180)*unitVect,'r','linewidth',2)
    % end
    figure()
    hlines=streamline(Xgrid,Ygrid,uField,vField,startx,starty);
    set(hlines,'LineWidth',2,'Color','r')
    hold on
    for j = 1:nairfoils

        plot(x(:,j),y(:,j),'-','color',rand(1,3),'linewidth',2)
        axis equal
        axis([-1 dist(end,1)+1.5 min(dist(:,2))-0.5 max(dist(:,2))+0.5])
        hold on
        xlabel('x/c')
        ylabel('y/c')
        quiver(xVel, yVel, cos(alpha(1)*pi/180)*unitVect, sin(alpha(1)*pi/180)*unitVect,'r','linewidth',2)
    end
    figure()
    quiver(Xgrid,Ygrid,uField,vField)
    hold on
    for j = 1:nairfoils

        plot(x(:,j),y(:,j),'-','color',rand(1,3),'linewidth',2)
        axis equal
        axis([-1 dist(end,1)+1.5 min(dist(:,2))-0.5 max(dist(:,2))+0.5])
        hold on
        xlabel('x/c')
        ylabel('y/c')
        quiver(xVel, yVel, cos(alpha(1)*pi/180)*unitVect, sin(alpha(1)*pi/180)*unitVect,'r','linewidth',2)
    end
    % figure()
    % UUU = sqrt(uField.^2 + vField.^2);
    % h=surf(Xgrid, Ygrid, UUU)
    % view(0,90)
    

end
>>>>>>> 2d703ea68bddf16c7d5ba0eb7850dd83c14bc9ee

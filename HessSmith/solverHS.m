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
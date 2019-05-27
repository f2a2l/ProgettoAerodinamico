function StartUI()

% Hess-Smith panel method for multi-airfoil case: the software allows to consider
% single or multi airfoil case, with an indefinite number of airfoils. For
% each airfoil, IGP parameters, the angle of attack, the leading edge
% coordinates (in a reference system defined by the body axes of the first airfoil
% and with the origin at its leading edge) and the relative chord must be specified (interactive
% input).
% Please remember: when inserting airfoil parameters, do include the square brackets (as in a vector definition).
% Original author: < Federico Messanelli: federico.messanelli@polimi.it >


  
AirfNumb = input([newline 'How many airfoils do you want to consider?' newline]);

% preallocation
alpha = zeros(1, AirfNumb);
aname = zeros(AirfNumb, 8);
if AirfNumb > 1
    dist = zeros(AirfNumb-1, 2);
    crel = zeros(1, AirfNumb-1);
end

disp([newline 'Please enter the (row) vector of parameters of the airfoil; do include square brackets.' newline])
for ii = 1:AirfNumb
    disp(['Insert profile ' int2str(ii)])
    disp('----------------')
    aname(ii,:) = input('Profile parameters: \n');
    alpha(ii) = input('Angle of attack: \n');
    if ii > 1
      dist(ii-1,1) = input('x position of leading edge: \n');
      dist(ii-1,2) = input('y position: \n');
      crel(ii-1) = input('Relative chord wrt first profile: \n');
    end
    disp(newline)
end



% -------------------------------------------------------------------------
% Hess-Smith method (external function)

if AirfNumb == 1
    [Cl, Cd, Cp, maxdCp, x, y, p, p1, SOL] = solverHS(80, aname, alpha, true);
else
    [Cl, Cd, Cp, maxdCp, x, y, p, p1, SOL] = solverHS(80, aname, alpha, dist, crel, true);
end



% -------------------------------------------------------------------------
% post processing and plot 
    
% nairfoils = AirfNumb;
% x0 = -1.1;
% yVel = [min(dist(:,2))-1:0.25:max(dist(:,2))+1]';
% xVel = -1.1*ones(length(yVel),1);
% unitVect = ones(length(yVel),1);
% for j = 1:nairfoils
%     
%     figure(99)
%     plot(x(:,j),y(:,j),'-','color',rand(1,3),'linewidth',2)
%     axis equal
%     axis([-1 max(dist(:,1))+1 min(dist(:,2))-1 max(dist(:,2))+1])
%     hold on
%     xlabel('x/c')
%     ylabel('y/c')
%     quiver(xVel, yVel, cos(alpha(1)*pi/180)*unitVect, sin(alpha(1)*pi/180)*unitVect,'r','linewidth',2)
% 
% end
%   
% 
% 

if AirfNumb == 1
    [uField,vField, Xgrid, Ygrid] = velocityField(p1, p1, alpha, 1, SOL, x, y);
else
    [uField,vField, Xgrid, Ygrid] = velocityField(p, p1, alpha, 1, SOL, x, y);
end

starty = (-1.5:0.025:1.5)';
startx = -0.48*ones(length(starty),1);
  
figure('Name', 'Flow field around multi element airfoil')
hlines=streamline(Xgrid,Ygrid,uField,vField,startx,starty);
set(hlines,'LineWidth',2,'Color','r')
hold on

strCl = [];
strCd = [];

for j = 1:AirfNumb

    plot(x(:,j),y(:,j),'-','color',rand(1,3),'linewidth',2)
    axis equal
    
    strCl = [strCl num2str(Cl(j))];
    strCd = [strCd num2str(Cd(j))];

    if j~=AirfNumb
        strCl = [strCl ' '];
        strCd = [strCd ' '];
    end

    title([newline 'CL: ' strCl newline 'CD: ' strCd newline newline])

end


% 
%   
% figure()
% quiver(Xgrid,Ygrid,uField,vField)
% hold on
% for j = 1:nairfoils
% 
%     plot(x(:,j),y(:,j),'-','color',rand(1,3),'linewidth',2)
%     axis equal
%     axis([-1 dist(end,1)+1.5 min(dist(:,2))-0.5 max(dist(:,2))+0.5])
%     hold on
%     xlabel('x/c')
%     ylabel('y/c')
%     quiver(xVel, yVel, cos(alpha(1)*pi/180)*unitVect, sin(alpha(1)*pi/180)*unitVect,'r','linewidth',2)
% 
% end

end
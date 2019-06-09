function StartHS()
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

disp(newline)
for ii = 1:AirfNumb
    disp(['Insert profile ' int2str(ii)])
    disp('----------------')
    aname(ii,:) = input('Please enter the (row) vector of parameters of the airfoil - do include square brackets: \n');
    alpha(ii) = input('Angle of attack: \n');
    if ii > 1
        disp('Please keep in mind that x and y axes are attached to (hence rotate with) the first profile.')
        dist(ii-1,1) = input('x position of leading edge wrt trailing edge of the first profile: \n');
        dist(ii-1,2) = input('y position: \n');
        crel(ii-1) = input('Relative chord wrt first profile: \n');
    end
    disp(newline)
end
disp('Do you want to show nodes on the airfoil surface on the plot?')
showNodes = input('Please type 1 for yes, 0 for no: ');
disp(newline)



% -------------------------------------------------------------------------
% Hess-Smith method (external function)

tic
if AirfNumb == 1
    [Cl, Cd, ~, ~, ~, x, y, p, p1, SOL] = solverHS(80, aname, alpha, true);
else
    [Cl, Cd, ~, ~, ~, x, y, p, p1, SOL, metaPan] = solverHS(80, aname, alpha, dist, crel, true);
end
toc


% -------------------------------------------------------------------------
% post processing and plot 


if AirfNumb == 1
    [uField,vField, Xgrid, Ygrid] = velocityField(p, AirfNumb, alpha, 1, SOL, x, y);
else
    [uField,vField, Xgrid, Ygrid] = velocityField(p, metaPan, AirfNumb, alpha, 1, SOL, x, y);
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

    plot(x{j},y{j},'-','color',rand(1,3),'linewidth',2)
    if showNodes
        scatter(x{j},y{j},'.','black')
    end
    axis equal
    
    strCl = [strCl num2str(Cl(j))];
    strCd = [strCd num2str(Cd(j))];

    if j~=AirfNumb
        strCl = [strCl '   '];
        strCd = [strCd '   '];
    end

    title([newline 'CL (left to right):   ' strCl newline 'CD (left to right):   ' strCd newline newline])

end

end
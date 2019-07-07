function plotStream(AirfNumb, p, x, y, metaPan, alpha, SOL, Cl, Cd, showNodes)

    if AirfNumb == 1
        [uField, vField, Xgrid, Ygrid] = velocityField(p, AirfNumb, alpha, 1, SOL, x, y);
    else
        [uField, vField, Xgrid, Ygrid] = velocityField(p, AirfNumb, alpha, 1, SOL, x, y, metaPan);
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
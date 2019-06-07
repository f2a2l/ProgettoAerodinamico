function [t] = lapTime_objFun(param)
    
    %% Unpacking param into Optimization Parameters
    c1_main = param(1,1);
    c2_main = param(1,2);
    c3_main = param(1,3);
    c4_main = param(1,4);
    x_t_main = param(1,5); 
    T_main = param(1,6); 
    rho_main = param(1,7); 
    beta_TE_main = param(1,8); 
    c1_flap = param(1,9); 
    c2_flap = param(1,10); 
    c3_flap = param(1,11); 
    c4_flap = param(1,12); 
    x_t_flap = param(1,13); 
    T_flap = param(1,14); 
    rho_flap = param(1,15); 
    beta_TE_flap = param(1,16); 
    AoA_main = param(1,17); 
    AoA_flap = param(1,18); 
    x_flap = param(1,19); 
    y_flap = param(1,20); 
    c_flap = param(1,21); 
    
    %% Generate Main airfoil and flap Geometry
    
    arflPar = [c1_main,c2_main,c3_main,c4_main,x_t_main,T_main,rho_main,beta_TE_main;...
				c1_flap,c2_flap,c3_flap,c4_flap,x_t_flap,T_flap,rho_flap,beta_TE_flap];
    
    
    %% Hess-Smith Calculation
    n_points = 80;
    
    [Cl,Cd,totLength,~,~,~,~,~,~,~] = solverHS(n_points,arflPar,[AoA_main AoA_flap],[x_flap y_flap],c_flap)

    %% XFOIL correction
    %TODO
    
    %% Check Quality of results
    if min(min(Cl,Cd)) > -1 && isreal(Cl) && isreal(Cd)
    
    %% 2D to 3D Correction
    b = 1.010; %wing span
    h = 0.67; %endplate height
    lambda = b / totLength;
    [Cl_new,Cd_new] = FWcorrections(Cl,Cd,lambda,b,h)
    fig=0;
    
    %% Lap Time Calculation
    [t] = sector(Cl_new, Cd_new,fig);
    
    
    else
        t = 10000;
    end
    
end
function [t] = lapTime_objFun(c1_main,c2_main,c3_main,c4_main,x_t_main,T_main,rho_main,beta_TE_main,...
                 c1_flap,c2_flap,c3_flap,c4_flap,x_t_flap,T_flap,rho_flap,beta_TE_flap,...
                 AoA_main, AoA_flap, x_flap, y_flap, c_flap)
    
    
    %% Generate Main airfoil and flap Geometry
    
    arflPar = [c1_main,c2_main,c3_main,c4_main,x_t_main,T_main,rho_main,beta_TE_main;
				c1_flap,c2_flap,c3_flap,c4_flap,x_t_flap,T_flap,rho_flap,beta_TE_flap];
            
    %% Hess-Smith Calculation
    n_points = 80;
    
    [Cl,Cd,totLength,~,~,~,~,~,~,~] = solverHS(n_points,arflPar,[AoA_main AoA_flap],[x_flap y_flap],c_flap);

    %% XFOIL correction
    %TODO
    
    %% 2D to 3D Correction
    b = 1.010; %wing span
    h = 0.67; %endplate height
    Lambda = b / totLength;
    [Cl_new,Cd_new] = FWcorrections(Cl,Cd,Lambda,b,h);
    
    %% Lap Time Calculation
    [t] = sector(Cl_new, Cd_new);
    
end
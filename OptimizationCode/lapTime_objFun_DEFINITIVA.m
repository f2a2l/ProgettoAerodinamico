function [t] = lapTime_objFun_DEFINITIVA(param)

% LAPTIME_OBJFUN Returns the time required for the Baku City Circuit 2nd Straight
%
% 1) It generates 2 airfoils reading the vecotr "param" which contains
% geometry info. #See documentation for arflPar generation and meaning.
%
% 2) Aerodynamic Performance is based on an inviscid preliminary estimation
% usign Hess-Smith solver, both for DRS-OFF and DRS-ON configuration.
%
% 3) A viscous correction is applied --> STILL TODO
%
% 4) At last, 2D to 3D geometry correction is applied
%
% 5) The code "sector.m" assembles the wing on a dynamic vehicle modeling to
% estimate acceleration, max speed and braking performances and returs the t.

tic;   
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
    
    %[Cl,Cd,totLength,~,~,~,~,~,~,~] = solverHS(n_points,arflPar,[AoA_main AoA_flap],[x_flap y_flap],c_flap);
    %[Cl_DRS,Cd_DRS,totLength,~,~,~,~,~,~,~] = solverHS_DRS(n_points,arflPar,[AoA_main AoA_flap],[x_flap y_flap],c_flap);
    %keyboard
    problem = solverVHS(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);
    problem_DRS = solverVHS_DRS(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);
    %keyboard
    %delta_cp = problem.maxdCp;
    cp = problem.Cp;
    cp_DRS = problem_DRS.Cp;
    
    %Stall check
    stallCheck;
    %[stall_or_not,~] = stall(delta_cp, [], c_flap);
    %keyboard
    if stall_flag == 0
%keyboard    
    plotCf = false;    
        [Cl_noDRS, Cd_noDRS, ~, ~] = problem.getBL(Re, plotCf);
        [Cl_DRS, Cd_DRS, ~, ~] = problem_DRS.getBL(Re,plotCf);
    %keyboard % Serve la parte con DRS
        Cl = [Cl_noDRS, Cl_DRS];
        Cd = [Cd_noDRS, Cd_DRS];
    %keyboard
    %% XFOIL correction
    %TODO
    %keyboard
    %% Check Quality of results
    %if min(min(min(Cl,Cd))) > -1 && isreal(Cl) && isreal(Cd)
        if Cd(1) > 0 && Cd(2) > 0 && isreal(Cl) && isreal(Cd)
    
    %% 2D to 3D Correction
        totLength = problem.xmax;
        b = 1.010; %wing span
        h = 0.67; %endplate height
        lambda = b / totLength;
        [Cl_new,Cd_new] = FWcorrections(Cl,Cd,lambda,b,h);
        fig=0;
    
    %% Lap Time Calculation
        [t] = sector(Cl_new, Cd_new,fig);
    
        %Try to sort out some problems. This is temporary. Find the error!
            if ~isreal(t)
                t = 10000;
            end
        
    %Penalization for nonphysical results and other problems
        else
            t = 10000;
        end
    
    % if stalled 
    else
        t = 10000;
    end
toc;    
end
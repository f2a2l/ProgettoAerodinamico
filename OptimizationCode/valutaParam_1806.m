% Min and Max values for the optimization procedure parameters
clear all
close all
clc
%% C1 - x position of 1st Mean Line control point
min_c1 = [0.25 0.3]; %default = 0.01;
max_c1 = [0.35 0.35]; %default = 0.96;
step_c1 = [0.1 0.1];

%% C2 - x position of 2nd Mean Line control point
min_c2 = [0.5 0.5]; %default = 0.02;
max_c2 = [0.8 0.8]; %default = 0.97;
step_c2 = [0.1 0.1];

%% C3 - y position of 1st Mean Line control point
min_c3 = [0 0.1]; %default = -0.074;
max_c3 = [0.3 0.3]; %default = 0.247;
step_c3 = [0.1 0.1];

%% C4 - y position of 2nd Mean Line control point
min_c4 = [0 0.15]; %default = -0.102;
max_c4 = [0.3 0.3]; %default = 0.206;
step_c4 = [0.1 0.1];

%% x_t - chordwise location of maximum thickness
min_x_t = [0.3 0.2]; %default = 0.2002;
max_x_t = [0.4 0.4]; %default = 0.4813;
step_x_t = [0.1 0.1];

%% T - maximum thickness
min_T = [0.1 0.1]; %default = 0.0246;
max_T = [0.18 0.15]; %default = 0.3227;
step_T = [0.1 0.1];

%% rho_bar - adimensional rho -> leading edge radius
min_rho = [0.5 0.175];
max_rho = [1 0.5]; %default = 1.4944;
step_rho = [0.1 0.1];

%% beta_TE_bar - trailing edge boat-tail angle
min_beta_TE = [2 2]; %default = 0.1452;
max_beta_TE = [3 3]; %default = 4.8724;
step_beta_TE = [0.1 0.1];

%% AoA main and flap - Angle of attack
min_AoA = [-5 30];
max_AoA = [25 70];
step_AoA = [1 -1];

%% Slot Definition -- CHECK AXES DEFINITION OF IGP
min_x = -0.02;%default = -0.02;
max_x = 0.02; % default = 0.016;
step_x = 0.01;

min_y = -0.015;
max_y = -0.03;
step_y = 0.01;

%% Flap Chord - (fraction of main airfoil chord, namely c = 1)
min_c_flap = 0.35;
max_c_flap = 0.5;
step_c_flap = 0.1;



%NOTE: If DRS on, AoA of the Main profile does not change. AoA of the flap
%tends towars zero-drag incidence. X distance is the same, Y distance
%becomes --> Y = 0.095 m.
duo = 1;
if duo
    
arflPar_min = [min_c1(1),min_c2(1),min_c3(1),min_c4(1),min_x_t(1),min_T(1),min_rho(1),min_beta_TE(1);...
				min_c1(2),min_c2(2),min_c3(2),min_c4(2),min_x_t(2),min_T(2),min_rho(2),min_beta_TE(2)];
[x_main, y_main, totLength] = multiGeometry(80, arflPar_min, [min_AoA(1) min_AoA(2)], [min_x min_y], min_c_flap,true);            
          
arflPar_max = [max_c1(1),max_c2(1),max_c3(1),max_c4(1),max_x_t(1),max_T(1),max_rho(1),max_beta_TE(1);...
				max_c1(2),max_c2(2),max_c3(2),max_c4(2),max_x_t(2),max_T(2),max_rho(2),max_beta_TE(2)];
           
[x_main, y_main, totLength] = multiGeometry(80, arflPar_max, [max_AoA(1) max_AoA(2)], [max_x max_y], max_c_flap,true);

else
arflPar = [0.3226,0.7107,0.2995,0.2961,0.2056,0.15,0.9343,2.9769;...
				0.2142,0.7558,0.2785,0.2652,0.3460,0.0901,0.9979,2.7961];

           
arflPar = [0.35,0.8,0.3,0.3,0.25,0.15,0.9343,2.9769;...
				0.2142,0.7558,0.2785,0.2652,0.3460,0.0901,0.9979,2.7961]; 
[x_main, y_main, totLength] = multiGeometry(80, arflPar, [max_AoA(1) max_AoA(2)], [max_x max_y], max_c_flap,true);           
end
% Min and Max values for the optimization procedure parameters

%% C1 - x position of 1st Mean Line control point
min_c1 = [0.1 0.1]; %default = 0.01;
max_c1 = [0.5 0.5]; %default = 0.96;
step_c1 = [0.1 0.1];

%% C2 - x position of 2nd Mean Line control point
min_c2 = [0.5 0.5]; %default = 0.02;
max_c2 = [0.8 0.8]; %default = 0.97;
step_c2 = [0.1 0.1];

%% C3 - y position of 1st Mean Line control point
min_c3 = [0 0]; %default = -0.074;
max_c3 = [0.3 0.3]; %default = 0.247;
step_c3 = [0.1 0.1];

%% C4 - y position of 2nd Mean Line control point
min_c4 = [0 0]; %default = -0.102;
max_c4 = [0.3 0.3]; %default = 0.206;
step_c4 = [0.1 0.1];

%% x_t - chordwise location of maximum thickness
min_x_t = [0.2 0.2]; %default = 0.2002;
max_x_t = [0.5 0.5]; %default = 0.4813;
step_x_t = [0.1 0.1];

%% T - maximum thickness
min_T = [0.01 0.01]; %default = 0.0246;
max_T = [0.15 0.15]; %default = 0.3227;
step_T = [0.1 0.1];

%% rho_bar - adimensional rho -> leading edge radius
min_rho = [0.175 0.175];
max_rho = [1 1]; %default = 1.4944;
step_rho = [0.1 0.1];

%% beta_TE_bar - trailing edge boat-tail angle
min_beta_TE = [1 1]; %default = 0.1452;
max_beta_TE = [3 3]; %default = 4.8724;
step_beta_TE = [0.1 0.1];

%% AoA main and flap - Angle of attack
min_AoA = [-5 -30];
max_AoA = [25; -70];
step_AoA = [1 -1];

%% Slot Definition -- CHECK AXES DEFINITION OF IGP
min_x = -0.01;%default = -0.02;
max_x = 0.02; % default = 0.016;
step_x = 0.01;

min_y = 0.01;
max_y = 0.016;
step_y = 0.1;

%% Flap Chord - (fraction of main airfoil chord, namely c = 1)
min_c_flap = 0.2;
max_c_flap = 0.4;
step_c_flap = 0.1;


%% ------------------------------------------------------------------ %%
%% ----------------- NON MODIFICARE OLTRE --------------------------- %%
%% ------------------------------------------------------------------ %%


%% IGP - Main Profile
c1_main = [min_c1(1):step_c1(1):max_c1(1)];
c2_main = [min_c2(1):step_c2(1):max_c2(1)];
c3_main = [min_c3(1):step_c3(1):max_c3(1)];
c4_main = [min_c4(1):step_c4(1):max_c4(1)];
x_t_main = [min_x_t(1):step_x_t(1):max_x_t(1)];
T_main = [min_T(1):step_T(1):max_T(1)];
rho_main = [min_rho(1):step_rho(1):max_rho(1)];
beta_TE_main = [min_beta_TE(1):step_beta_TE(1):max_beta_TE(1)];

%% IGP - Flap Profile
c1_flap = [min_c1(2):step_c1(2):max_c1(2)];
c2_flap = [min_c2(2):step_c2(2):max_c2(2)];
c3_flap = [min_c3(2):step_c3(2):max_c3(2)];
c4_flap = [min_c4(2):step_c4(2):max_c4(2)];
x_t_flap = [min_x_t(2):step_x_t(2):max_x_t(2)];
T_flap = [min_T(2):step_T(2):max_T(2)];
rho_flap = [min_rho(2):step_rho(2):max_rho(2)];
beta_TE_flap = [min_beta_TE(2):step_beta_TE(2):max_beta_TE(2)];

%% AoA
AoA_main = [min_AoA(1):step_AoA(1):max_AoA(1)];
AoA_flap = [min_AoA(2):step_AoA(2):max_AoA(2)];

%% Slot Definition
x = [min_x:step_x:max_x];
y = [min_y:step_y:max_y];

%% Flap Chord (relative to Main Profile Chord c = 1)
c_flap = [min_c_flap:step_c_flap:max_c_flap];


%NOTE: If DRS on, AoA of the Main profile does not change. AoA of the flap
%tends towars zero-drag incidence. X distance is the same, Y distance
%becomes --> Y = 0.095 m.


%% To optimization Solver
lb = [min_c1(1),min_c2(1),min_c3(1),min_c4(1),min_x_t(1),min_T(1),min_rho(1),min_beta_TE(1),...
       min_c1(2),min_c2(2),min_c3(2),min_c4(2),min_x_t(2),min_T(2),min_rho(2),min_beta_TE(2),... 
       min_AoA(1),min_AoA(2),...
       min_x,min_y,min_c_flap];
ub = [max_c1(1),max_c2(1),max_c3(1),max_c4(1),max_x_t(1),max_T(1),max_rho(1),max_beta_TE(1),...
       max_c1(2),max_c2(2),max_c3(2),max_c4(2),max_x_t(2),max_T(2),max_rho(2),max_beta_TE(2),... 
       max_AoA(1),max_AoA(2),...
       max_x,max_y,max_c_flap];
  
clearvars -except lb ub
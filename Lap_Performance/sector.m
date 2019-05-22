function [T_sector] = sector(CL, CD,fig)
% SECTOR is a function which takes ad inputs the
% aerodynamic coefficients of REAR WING --> from optimization loop
%
%         [T_sector] = sector(CL, CD, fig)
%
% The input structure is the following:
%       DRS OFF |  DRS ON 
% CL = [ CL(1)  |   CL(2) ]
% CD = [ CD(1)  |   CD(2) ]
%
% The output is the time to run the sector: T_sector
% Set fig == true for plot


% Generic inputs
model

% Specific inputs --> FIXED for Baku Circuit
L = 850;                 % [m]       sector length
D = 160;                 % [m]       distance where DRS can be activated
u_in = 95 / 3.6;         % [m/s]     speed at turn 2
u_out = 100 / 3.6;       % [m/s]     speed at turn 3

disp('**********************************************************')
disp('* BAKU CITY CIRCUIT - straight between turn 2 and turn 3 *')
disp('**********************************************************')
disp(['Straight length:     ', num2str(L), ' m'])
disp(['DRS activatio point: ', num2str(D), ' m'])
disp(['Entry speed:         ', num2str(u_in*3.6), ' km/h'])
disp(['Exit speed:          ', num2str(u_out*3.6), ' km/h'])


% Aerodynamic coefficient of F1 CAR
CL = CL_car + CL;                
CD = CD_car + CD; 


% Compute maximum speed
u_max = max_speed(CD(2),fig);   % [m/s]    maximum car speed

fprintf('\n')
disp('-----------------------------------------------')
disp(['| Maximum speed: ',num2str(u_max),' m/s or ',...
    num2str(u_max*3.6),' km/h |'])
disp('-----------------------------------------------')
fprintf('\n')

%% Sub sectors

%%% Acceleration - DRS OFF %%%
% Compute time to reach DRS activation point at a distance D from turn 2
% and compute speed when DRS is activated
[t_SB1,u_SB1] = acceleration_DRSOFF(u_in, D, CL(1), CD(1),fig);
s_SB1 = D;

%%% Acceleration - DRS ON %%%
% Compute space and time needed to reach max speed 
% from DRS activation point
[s_SB2, t_SB2] = acceleration_DRSON(u_SB1, u_max, CL(2),CD(2),fig);

%%% Breaking - DRS OFF %%%
% Compute first time and distance needed to perform a braking from
% maximum car speed to turn 3 speed
[s_SB4, t_SB4] = braking(u_max,u_out,CL(1),CD(1),fig);

%%% Max speed - DRS ON %%%
% First compute the track length left then the time needed to run this
% part maximum speed
s_SB3 = L - D - s_SB2 - s_SB4;
t_SB3 = s_SB3 / u_max;

% Sector time is the sum of all sub sector times
%    SB1    SB2      SB3     SB4
T_sector = t_SB1 + t_SB2 + t_SB3 + t_SB4;

%% Display results
fprintf('\n')
disp('****** SUB SECTOR 1 ******')
disp( ['Acceleration time: ',num2str(t_SB1,4),' s'])
disp( ['Acceleration final speed (NO DRS): ',num2str(u_SB1,4),' m/s'])
fprintf('\n')
disp('****** SUB SECTOR 2 ******')
disp( ['Acceleration time: ',num2str(t_SB2),' s'])
disp( ['Acceleration distance: ',num2str(s_SB2),' m'])
fprintf('\n')
disp('****** SUB SECTOR 3 ******')
disp(['Max speed straigth distance: ', num2str(s_SB3),' m'])
disp(['Max speed straigth time: ', num2str(t_SB3),' s'])
fprintf('\n')
disp('****** SUB SECTOR 4 ******')
disp( ['Braking time: ',num2str(t_SB4),' s'])
disp( ['Braking distance: ',num2str(s_SB4),' m'])
fprintf('\n')

disp('--------------------------')
disp(['| Sector time: ',num2str(T_sector),' s |']) 
disp('--------------------------')

end
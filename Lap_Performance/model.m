%% Model

% This scipt contains the data used to compute sector performance

g   =  9.81;          % [m/s^2] gravitational field intensity
rho = 1.225;          % [kg/m^3] air density at sea level 
mu  =   1.5;          % [-] grip level (maximum value: 1.6)
m   =   823;          % [kg] total mass of F1 car
Sa   =   1.3;         % [m^2] typical frontal area

% Typical values of Sc --> surface * aerodynamic coefficient

Scx = 1.2;            % [m^2] longitudinal direction 
CD_car = .4*(Scx/Sa); % --> 40% of car contribution without rear wing

Scz =   6;            % [m^2] vertical direction
CL_car = .4*(Scz/Sa); % --> 40% of car contribution without rear wing
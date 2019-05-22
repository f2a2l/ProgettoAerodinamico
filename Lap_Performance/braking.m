function [s_finale, t_finale] = braking(u_in,u_out,CL,CD,fig)

% BRAKING is a function which allows to compute space and time
% needed for braking event, given the initial speed, the final speed and
% lift/drag coefficients of the rear wing:
%
%           [s_finale, t_finale] = braking(u_in,u_out,CL,CD,fig)
%
% Inputs:
% - u_in       speed at which braking event begins
% - u_out      turn speed --> constant along the turn
% - CL         lift coefficient of rear wing (DRS OFF)
% - CD         drag coefficient of rear wing (DRS OFF)
%
% Set fig == true for plot

% Import data for braking
model

% Aerodynamic parameters
zeta = rho * Sa * CL / 2; 
eta   = rho * Sa * CD / 2; 

% Initialize variables
t_step = 0;
s = 0;
u = u_in;

T = 0;
S = 0;
U = u_in;
dt = 0.0001;

% Define maximum braking deceleration
f_max =@(u) mu * (g + (zeta)/m * u.^2 ) + eta / m * u.^2;

% Imposition of maximum deceleration (g) 
F_lim = 5 * 9.81; % [m/s^2]
F_plot = F_lim;

while (u >= u_out)
   
    t_step = t_step + dt;
    F = f_max(u);
    
    if F >= F_lim
        F = F_lim;
    end
    u_new = -F * t_step + u;
    s_new = -0.5*F * t_step^2 + u * t_step + s; 
    T = [T, t_step];
    U = [U, u_new];
    S = [S, s_new];
    u = u_new;
    s = s_new;
    F_plot=[F_plot,F];

end


t_finale = T(end) * length(T);
s_finale = S(end);

% disp( ['Braking time: ',num2str(t_finale),' s'])
% disp( ['Braking distance: ',num2str(s_finale),' m'])

%% Plot

if fig == true
figure()
subplot(3,1,1)
plot(linspace(0,t_finale,length(U)),U)
title (['SUB SECTOR 4 - Braking Performance from ',num2str(u_in,3),' m/s to '...
    ,num2str(u_out,3), ' m/s'])
ylabel 'Speed [m/s]'
xlabel 'Time [s]'
hold on 
grid on
plot(0, u_in,'o','LineWidth',1)
plot(t_finale,u_out,'o','LineWidth',1)
legend('u(t)','u_{max}','u_{turn}','Location','best')

subplot(3,1,2)
plot(linspace(0,t_finale,length(S)),S)
ylabel 'Distance [m]'
xlabel 'Time [s]'
grid on
hold on
plot(0,0,'o','LineWidth',1)
plot(t_finale,S(end),'o','LineWidth',1)

subplot(3,1,3)
plot(linspace(0,t_finale,length(F_plot)),F_plot./9.81)
ylabel 'Deceleration'
xlabel 'Time [s]'
ylim([0, 6])
grid on
hold on
plot(0,F_plot(1),'o','LineWidth',1)
plot(t_finale,F_plot(end),'o','LineWidth',1)

end
end
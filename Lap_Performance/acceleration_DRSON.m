function [s_finale, t_finale] = acceleration_DRSON(u_in, u_max,CL,CD,fig)

% ACCELERATION_DRSON is a function which allows to compute space and time 
% needed for acceleration event when the DRS is open, from an initial
% speed u_in to the maximum car speed u_max
%
%           [s_finale, t_finale] = acceleration_DRSON(u_in, u_max,CL,CD,fig)
%
% Inputs:
% - u_in       speed at which acceleration with DRS ON begins
% - u_max      maximum F1 car speed
% - CL         lift coefficient of rear wing (DRS ON)
% - CD         drag coefficient of rear wing (DRS ON)
%
% Set fig == true for plot

% Import data for braking
model

% Aerodynamic parameters
zeta = rho * Sa * CL / 2; 
eta   = rho * Sa * CD / 2; 

% Typical F1 values
% 0 - 100 kmh --> 2.4 s
% 0 - 200 kmh --> 4.4 s
% 0 - 300 kmh --> 8.4 s

t = 0;
s = 0;
u = u_in;

T = 0;
S = 0;
U = u_in;
dt = 0.0001;

% Define maximum acceleration
acc_max =@(u) mu * (g + (zeta)/m * u.^2 ) + eta / m * u.^2;

% Imposition of maximum acceleration (g)
a_lim = 1.1 * 9.81; % [m/s^2]
A_plot = acc_max(u_in);

while (u <= u_max)
   
    t = t + dt;
    A = acc_max(u);
    
    if A >= a_lim
        A = a_lim;
    end
    u_new = A * dt + u;
    s_new = 0.5*A * dt^2 + u * dt + s; 
    T = [T, t];
    U = [U, u_new];
    S = [S, s_new];
    u = u_new;
    s = s_new;
    A_plot = [A_plot,A];

end

t_finale = T(end);
s_finale = S(end);

% disp( ['Acceleration time: ',num2str(t_finale),' s'])
% disp( ['Acceleration distance: ',num2str(s_finale),' m'])

%% Plot

if fig == true
figure()
subplot(3,1,1)
plot(linspace(0,t_finale,length(U)),U)
title (['SUB SECTOR 2 - Acceleration Performance from ',num2str(u_in,3),' m/s to '...
    ,num2str(u_max,3), ' m/s'])
ylabel 'Speed [m/s]'
xlabel 'Time [s]'
hold on 
grid on
plot(0, u_in,'o','LineWidth',1)
plot(t_finale,u_max,'o','LineWidth',1)
legend('u(t)','u_{turn}','u_{max}','Location','best')

subplot(3,1,2)
plot(linspace(0,t_finale,length(S)),S)
ylabel 'Distance [m]'
xlabel 'Time [s]'
grid on
hold on
plot(0,0,'o','LineWidth',1)
plot(t_finale,S(end),'o','LineWidth',1)

subplot(3,1,3)
plot(linspace(0,t_finale,length(A_plot)),A_plot./9.81)
ylabel 'Acceleration'
xlabel 'Time [s]'
ylim([0, 3])
grid on
hold on
plot(0,A_plot(1),'o','LineWidth',1)
plot(t_finale,A_plot(end),'o','LineWidth',1)
end

end 
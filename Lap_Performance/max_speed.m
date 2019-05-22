function [u_max]=max_speed(CD)

% MAX_SPEED is function which allows to compute maximum speed of an
% F1 car, given the car drag coefficient. The maximum engine power
% is considered fixed for the Ferrari SF90: 990 CV.
%
%           u_max=max_speed(CD)


% source: WIKIPEDIA Ferrari SF90
max_CV = 990;                  % [Hp] maximum horse power
max_power = max_CV/1.36*1000;  % [kW]

% General inputs
Sa = 1.3; % [m^2]
rho = 1.225; 

u_max = (2*max_power / (rho*Sa*CD) )^(1/3);

D_power =@(V) 0.5 * rho * Sa * CD * V.^3; % [kW]
V = linspace(0,110,100);

figure()
subplot(2,1,1)
plot(V,D_power(V)./1000,'LineWidth',1)
grid on
hold on
plot(V,max_power*ones(1,length(V))./1000,'LineWidth',1)
hold on
plot(u_max,D_power(u_max)./1000,'ob','LineWidth',1.5)
title(['Maximum speed: ',num2str(u_max,3), ' m/s or ',...
    num2str(u_max*3.6,3), ' km/h'])
xlabel('Speed [m/s]')
ylabel('Power [kW]')
xlim([0, V(end)])
legend('Drag Power','Max Engine Power','Location','NorthWest')

subplot(2,1,2)
plot(V.*3.6,D_power(V)./1000,'LineWidth',1)
grid on
hold on
plot(V.*3.6,max_power*ones(1,length(V))./1000,'LineWidth',1)
hold on
plot(u_max.*3.6,D_power(u_max)./1000,'ob','LineWidth',1.5)
xlabel('Speed [km/h]')
ylabel('Power [kW]')
xlim([0, V(end)*3.6])
legend('Drag Power','Max Engine Power','Location','NorthWest')


end
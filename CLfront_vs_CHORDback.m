close all
clear
clc

% Front airfoil geometry parameters
arflPar_single = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];

% Tandem configuration geometry parameters
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5;
           0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];

alpha_f = 5;                  % fixed AoA front airfoil
alpha_b = 30 : 5 : 60;        % range of AoA back airfoil
chord_b = 0.25 : 0.05 : 0.35; % range of chord back airfoil


Cl_double = zeros(2,1);
Cl_plot= zeros(length(chord_b),length(alpha_b)); % front airfoil CL matrix

[Cl, ~, ~, ~] = solverHS(80, arflPar_single, alpha_f); % reference value of CL front

for j = 1 : length(chord_b)
    for i = 1 : length(alpha_b)
        
        [Cl_double, ~, ~, ~] = solverHS(80, arflPar, [alpha_f, alpha_b(i)], [1.05, -0.05], chord_b(j));
        Cl_plot(j,i) = Cl_double(1);
        
    end 
end

%% Plot 

figure()
for j = 1 : length(chord_b)
    plot(alpha_b,Cl_plot(j,:),'-o','DisplayName',['chord_{back} = ',num2str(chord_b(j))])
    hold on
    grid on
end
legend('show','Location','NW')
title ('Variazione del CL_{front} al variare dell'' \alpha_{back}')
xlabel ('\alpha_{back}')
ylabel ('CL_{front}')
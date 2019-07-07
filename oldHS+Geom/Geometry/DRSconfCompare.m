% This is to compare DRS-OFF and DRS-ON configurations

clear all
clc
close all

set(0,'defaultAxesFontSize',18)
set(0,'defaultAxesTickLabelInterpreter','latex')
%set(0,'defaultLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')

%% Select IGP parameters to generate airfoils
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5;
					 0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
                 
                 
                 
%% Generate DRS-OFF Geometry                 
[x_off, y_off, ~] = multiGeometry(80, arflPar, [3 -30], [+0.05 +0.05], 0.3);

%% Generate DRS-ON Geometry
[x_on, y_on, ~] = multiGeometry_DRS(80, arflPar, [3 -30], [+0.05 +0.05], 0.3);


%% Plot Geometries compared

figure(1)
subplot(1,2,1)
plot(x_off,y_off)
%ylim([-0.3,0.4])
%xlim([0,1.305])
title('DRS-OFF')
axis equal

%figure(2)
subplot(1,2,2)
plot(x_on,y_on)
%ylim([-0.3,0.4])
%xlim([0,1.305])
title('DRS-ON')
axis equal
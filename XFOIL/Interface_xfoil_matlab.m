
clc
clear all

% Inputs: 

coord = 'NACA 0012';
alpha =  2;
Re    =  1e6; % Visc 'on' is automatically started if Re>0.
Mach  =  0;
Npanels = '300'; %'panels n 300';
Niterations = '50'; %'iter 50';

% Calling the fuction:
[polar,foil] = xfoil(coord,alpha,Re,Mach,Npanels, Niterations);

%Plot:
figure; subplot(3,1,[1 2]);
plot(foil(1).xcp(:,end),foil(1).cp(:,end)); 
xlabel('x');
ylabel('C_p'); 
title(sprintf('%s @ %g\\circ',polar.name,foil(1).alpha(end))); 
set(gca,'ydir','reverse');
subplot(3,1,3);
I = (foil(1).x(:,end)<=1); 
plot(foil(1).x(I,end),foil(1).y(I,end)); xlabel('x');
ylabel('y'); axis('equal'); 

% Output Data:
AoA = polar.alpha
Re = Re
CL = polar.CL
CD = polar.CD
CDp = polar.CDp %Parasite Drag (Skin Friction)
Cm = polar.Cm
Top_xtr = polar.Top_xtr % Location of transition flow on top.
Bot_xtr = polar.Bot_xtr % Location of transition flow on bottom.



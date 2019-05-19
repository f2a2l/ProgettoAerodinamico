% Plot profile
clc
clear all
close all

%% Profile Selection
fprintf(['Supported Profiles: \n\n' ...
                '* NACA 4 digits \n' ...
                '* NACA 5 digits \n' ...
                '* Benzing Profiles \n\n']);
aname = input('Profile: ','s');

if isempty(aname)
    aname = 'Naca0012';
end

%% Spacing Selection
fprintf(['\nSpacing Law: \n\n'...
                '[1] constant \n' ...
                '[2] halfcos \n' ...
                '[3] cos           <-- default \n\n']);
spacing = input('Spacing: ','s');

if isempty(spacing)
    spacing = 'cos';
end

%% Generate Profile Points
[x,y] = AirfoilShape(aname,spacing);
 

%% Plot

%Plot Settings
x0=10; %x position
y0=10; %y position
width=1800; %x size
height=300; %y size
set(gcf,'position',[x0,y0,width,height])



plot(x,y)
xlim([0,1])
%ylim([-0.2, 0.2])
text(0.3,mean(y),aname,'fontsize',20);
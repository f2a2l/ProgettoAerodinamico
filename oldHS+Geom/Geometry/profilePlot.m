%% Profile Selection
fprintf(['Please enter a set of parameters for the profile.' newline]);
aname = input('Vector of IGP parameters (please include square brackets): ');

notValidInput = false;
lng = 8;

try
    lng = length(aname);
catch
    notValidInput = true;
end

if isempty(aname)
    warning('Using default profile.')
    aname = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
elseif lng ~= 8 || notValidInput
    warning('Invalid input: switching to default')
    aname = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
end

np = input('Select number of points: ');

if isempty(np)
    warning('Using default number of points (80).')
    np = 80;
end

%% Generate Profile Points
[x,y] = AirfoilShape(aname, np);
 

%% Plot

figure('Name','Profile geometry')
% x0=10; %x position
% y0=10; %y position
% width=1800; %x size
% height=300; %y size
% set(gcf,'position',[x0,y0,width,height])

scatter(x,y, 'filled', 'black')
axis equal

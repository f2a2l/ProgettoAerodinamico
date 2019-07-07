% Creates a .dat file containing the points of the airfoils


%% Generate airfoil geometry from optimization results
if exitflag == 1
    
%Unpack optimalSolution vector
unpackOptimalSolution;

arflPar = [c1_main, c2_main, c3_main, c4_main, x_te_main, T_main, rho_main, beta_te_main;
            c1_flap, c2_flap, c3_flap, c4_flap, x_te_flap, T_flap, rho_flap, beta_te_flap];

[x_main, y_main, xtotLength,ytotLength] = multiGeometry(80, arflPar, [AoA_main AoA_flap], [x_flap y_flap], c_flap,true);
[x_DRS, y_DRS, ~, ytotLength_DRS] = multiGeometry_DRS(80, arflPar, [AoA_main AoA_flap], [x_flap y_flap], c_flap,true);

else

disp('Optimization Procedure FAILED!')
    
% Lines 20 - 26 are only meant to execute the code.
% For the final simulation, they have to be removed
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5;
					 0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];

[x_main, y_main, totLength] = multiGeometry(80, arflPar, [3 10], [-0.05 -0.05], 0.3, true);
[x_DRS, y_DRS, totLength] = multiGeometry_DRS(80, arflPar, [3 10], [-0.05 -0.05], 0.3, true);

end

%% Extract 2 vectors from arrays
x_flap = [x_main{2}];
y_flap = [y_main{2}];

x_main = [x_main{1}];
y_main = [y_main{1}];

x_DRS_main = [x_DRS{1}];
y_DRS_main = [y_DRS{1}];

x_DRS_flap = [x_DRS{2}];
y_DRS_flap = [y_DRS{2}];

%% Scaling airfoils to match Rear Wing geometry regulations

%Imposed Length
maxLength = 350;

%Scale factor
scaleFactor = 350 / xtotLength;

%Scaling DRS-OFF and DRS-ON conf. points
x_flap = x_flap * scaleFactor;
y_flap = y_flap * scaleFactor;

x_main = x_main * scaleFactor;
y_main = y_main * scaleFactor;

x_DRS_main = x_DRS_main * scaleFactor;
y_DRS_main = y_DRS_main * scaleFactor;

x_DRS_flap = x_DRS_flap * scaleFactor;
y_DRS_flap = y_DRS_flap * scaleFactor;

%% n-Points (necessary for the file header)
%n = length(x_main) + length(x_flap);
%npoints = num2str(n);

%% z-coordinate (set to 0, since 2D)
%z = zeros(n,1);

%% Generate .dat files
writeGeomFile(x_main,y_main,[],'main')
writeGeomFile(x_flap,y_flap,[],'flap')
writeGeomFile(x_DRS_flap,y_DRS_flap,[],'flap_DRS')
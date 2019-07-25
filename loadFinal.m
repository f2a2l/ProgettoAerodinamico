addPaths
load('opti_results.mat')
load('main_mod.mat')
unpackOptimalSolution
saveProfiles

%% best fit of Xfoil-corrected profile
% mp = findBestFitProfile(mainmod(:,1), mainmod(:,2), arflPar(1,:));
% arflPar(1,:) = mp;
% save('final_geom', 'arflPar')

%% load final geometry
load('final_geom.mat')
profile = solverVHS(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);
profile_DRS = solverVHS_DRS(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);


% set(0,'defaultAxesFontSize',18)
% set(0,'defaultAxesTickLabelInterpreter','latex')
% set(0,'defaultTextInterpreter','latex')


%% DRS OFF

% [x, y] = multiGeometry(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);
% 
% figure
% hold on
% for ii = 1:length(x)
%     xp = x{ii};
%     yp = y{ii};
%     plot(xp, yp, 'black', 'LineWidth', 2);
% end
% axis equal
% set(get(gca, 'XLabel'), 'String', 'x/c');
% set(get(gca, 'YLabel'), 'String', 'y/c');



%% DRS ON

% [x, y] = multiGeometry_DRS(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);
% 
% figure
% hold on
% for ii = 1:length(x)
%     xp = x{ii};
%     yp = y{ii};
%     plot(xp, yp, 'black', 'LineWidth', 2);
% end
% axis equal
% set(get(gca, 'XLabel'), 'String', 'x/c');
% set(get(gca, 'YLabel'), 'String', 'y/c');
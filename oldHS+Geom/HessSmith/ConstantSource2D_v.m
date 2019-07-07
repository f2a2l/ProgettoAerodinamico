function [u] = ConstantSource2D_v(sigma_J, panel_J, panels)


%R = Rotation (panel_J.beta);

x1 = -panel_J.d/2;
x2 = panel_J.d/2;

CPan = [panels.C];

P_I = panel_J.R*[[CPan.x]-panel_J.C.x; [CPan.y]-panel_J.C.y];

[upx, upy] = ConstantSource2D_local (sigma_J, x1, x2, P_I(1,:), P_I(2,:));

u = panel_J.R'*[upx; upy];

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [u, v] = ConstantSource2D_local (sigma, x1, x2, xp, yp)

u = (sigma/(4*pi)) * log (((xp-x1).^2 + yp.^2)./((xp-x2).^2 + yp.^2));
v = (sigma/(2*pi)) * (atan((xp-x1)./yp) - atan((xp-x2)./yp));

return


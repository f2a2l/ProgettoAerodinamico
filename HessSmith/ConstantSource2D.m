function [u] = ConstantSource2D (sigma_J, panel_J, panel_I)


%R = Rotation (panel_J.beta);

x1 = -panel_J.d/2;
x2 = panel_J.d/2;

P_I = panel_J.R*[panel_I.C.x-panel_J.C.x; panel_I.C.y-panel_J.C.y];

[upx, upy] = ConstantSource2D_local (sigma_J, x1, x2, P_I(1), P_I(2));

u = panel_J.R'*[upx; upy];

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [u, v] = ConstantSource2D_local (sigma, x1, x2, xp, yp)

u = (sigma/(4*pi)) * log (((xp-x1)^2 + yp^2)/((xp-x2)^2 + yp^2));
v = (sigma/(2*pi)) * (atan((xp-x1)/yp) - atan((xp-x2)/yp));

return


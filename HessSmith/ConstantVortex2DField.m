function [u, v] = ConstantVortex2DField (gamma_J, panel_J, x, y)

%

%R = Rotation (panel_J.beta);

x1 = -panel_J.d/2;
x2 = panel_J.d/2;

P_I = panel_J.R*[x'-panel_J.C.x; y'-panel_J.C.y];

[upx, upy] = ConstantVortex2D_local (gamma_J, x1, x2, P_I(1,:), P_I(2,:));

UU = panel_J.R'*[upx; upy];

u = UU(1,:);
v = UU(2,:);

return

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [u, v] = ConstantVortex2D_local (gamma, x1, x2, xp, yp)


u = (gamma/(2*pi)) * (atan((xp-x1)./yp) - atan((xp-x2)./yp));
v = -(gamma/(4*pi)) * log (((xp-x1).^2 + yp.^2)./((xp-x2).^2 + yp.^2));

return


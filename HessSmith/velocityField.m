function [u,v, Xgrid, Ygrid] = velocityField(p, p1, alpha, U, SOL, xa ,ya)


ntot = length(p.panel);
nairfoils = length(p1);
nterz = (ntot - nairfoils)/nairfoils;
alpha1 = alpha(1); % il sist. di riferimento � solidale al primo profilo
uInf1 = U*[cos(deg2rad(alpha1)); sin(deg2rad(alpha1))];

for i = 1:ntot
    xC(i) = p.panel(i).C.x;
    yC(i) = p.panel(i).C.y;
end



xmin = min(xC) - 0.5;
xmax = max(xC) + 0.5;
ymin = min(yC) - 0.5;
ymax = max(yC) + 0.5;

Ngrid = 150;

xgrid = linspace(xmin,xmax,Ngrid)';
ygrid = linspace(ymin,ymax,Ngrid)';



[Xgrid, Ygrid] = meshgrid(xgrid,ygrid);
xpoints = Xgrid(:);
ypoints = Ygrid(:);

u = zeros(1,Ngrid^2);
v = zeros(1,Ngrid^2);
% u = zeros(Ngrid, Ngrid);
% v = zeros(Ngrid, Ngrid);

for j = 1:ntot-nairfoils
  [us, vs] = ConstantSource2DField (SOL(j), p.panel(j), xpoints, ypoints);
  u = u + us ; 
  v = v + vs ;
end

for k = 1:nairfoils
  for j =  ((k-1)*nterz + 1 ):k*nterz
    [uv vv] = ConstantVortex2DField (SOL(end-nairfoils+k), p.panel(j), xpoints, ypoints);
    u = u + uv ; %+ ti*uv;
    v = v + vv ;
  end
end

u = u + uInf1(1);
v = v + uInf1(2);

for k = 1:nairfoils
  [in,on] = inpolygon(xpoints,ypoints,xa(:,k),ya(:,k));
  u(in | on) = 0;
  v(in | on) = 0;
end

u = reshape(u,Ngrid,Ngrid);
v = reshape(v,Ngrid,Ngrid);



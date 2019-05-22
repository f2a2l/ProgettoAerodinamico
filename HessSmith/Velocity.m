function[v] = Velocity (p, alpha, U, SOL)


n = length(p.panel);

uInf = U*[cos(deg2rad(alpha)); sin(deg2rad(alpha))];

v = zeros (n-1, 1);


ti = [cos([p.panel(1:n-1).beta]); sin([p.panel(1:n-1).beta])];

% preallocation
us = zeros(2,n-1);
uv = us;

for j = 1:n-1

  idx = 1:n-1;
  idx(j) = [];
  us(:,idx) = ConstantSource2D_v(SOL(j), p.panel(j), p.panel(idx));
  uv(:,idx) = ConstantVortex2D_v(SOL(end), p.panel(j), p.panel(idx));

  us(:,j) = p.panel(j).R'*[0; 1/2]*SOL(j);
  uv(:,j) = p.panel(j).R'*[1/2; 0]*SOL(end);
  v = v + sum(ti.*us,1)' + sum(ti.*uv,1)';

end

v = v + ti'*uInf;

return
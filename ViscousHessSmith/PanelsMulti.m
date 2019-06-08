function [p, metaPan] = PanelsMulti(p1)


nairfoil = length(p1);

p = struct;

% get no pts for each airfoil
nptz = zeros(1,nairfoil); 
for ii = 1:nairfoil
    nptz(ii) = length(p1(ii).panel);
end

% get number of TRUE panels (excludes Kutta panels) for each airfoil
npan = nptz - 1;

% zero index for panels of each airfoil
idx_zeroPan = zeros(1, nairfoil);
for ii = 2:nairfoil
    idx_zeroPan(ii) = sum(npan(1:(ii-1)));
end

% zero index for kutta panels
idx_zeroKutta = sum(npan);

% ntot = sum(nptzsa);

% nterz = (ntot-nairfoil)/nairfoil;

               
for k = 1:nairfoil

    % nterz = npoints - 1;

    for j = 1:npan(k)
        i = idx_zeroPan(k) + j;
        p.panel(i).P1.x  = p1(k).panel(j).P1.x;
        p.panel(i).P2.x  = p1(k).panel(j).P2.x;
        p.panel(i).P1.y  =  p1(k).panel(j).P1.y;
        p.panel(i).P2.y  = p1(k).panel(j).P2.y;
        p.panel(i).C.x   = p1(k).panel(j).C.x;
        p.panel(i).C.y   =  p1(k).panel(j).C.y;
        p.panel(i).d     = p1(k).panel(j).d;
        p.panel(i).beta  = p1(k).panel(j).beta ;
        p.panel(i).R  = p1(k).panel(j).R ;
    end

    i = idx_zeroKutta + k;
    npoints = nptz(ii);
    p.panel(i).P1.x  = p1(k).panel(npoints).P1.x;
    p.panel(i).P2.x  = p1(k).panel(npoints).P2.x;
    p.panel(i).P1.y  =  p1(k).panel(npoints).P1.y;
    p.panel(i).P2.y  = p1(k).panel(npoints).P2.y;
    p.panel(i).C.x   = p1(k).panel(npoints).C.x;
    p.panel(i).C.y   =  p1(k).panel(npoints).C.y;
    p.panel(i).d     = p1(k).panel(npoints).d;
    p.panel(i).beta  = p1(k).panel(npoints).beta ;
    p.panel(i).R  = p1(k).panel(npoints).R ;
end

metaPan = struct;
metaPan.npan = npan;
metaPan.idx_zeroPan = idx_zeroPan;

    
return

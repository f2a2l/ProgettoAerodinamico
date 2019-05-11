function [p] = PanelsMulti(p1)


nairfoil = length(p1);
npoints = length(p1(1).panel);

p = struct;

ntot = nairfoil*npoints;

nterz = (ntot-nairfoil)/nairfoil;

               
for k = 1:nairfoil

    for j = 1:(npoints-1)
        i = ((k-1)*nterz + j);
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

    i = ((nairfoil-1)*nterz + npoints -1) + k;
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


    
return

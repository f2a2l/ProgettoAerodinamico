function writeGeomFile(x,y,z,n,name)

% WRITEGEOMFILE writes a .dat file containing (x,y,z) coordinates
%
% This is meant to create a ready-to-use .dat file for CFD meshing in Pointwise
%
% Example:
%
%   writeGeomFile(x,y,z,20,naca0012)

name = strcat(name,'.dat');
npoints = num2str(n);

fid = fopen(name,'w');
fprintf(fid, npoints);

for k = 1:n
    
    fprintf(fid, '\n%f   %f   %f',x(k), y(k), z(k));
    
end

fclose(fid);

end
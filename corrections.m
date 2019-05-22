function[Cl_new,Cd_new] = corrections(Cl,Cd_0,Lambda,b,h,s,t)
% Cl = Cl before corrections
% Cd_0 = viscous drag
% Lambda = wing aspect ratio
% b = wing spanwise
% h = endplate's high
% s = parameter taken from pag. 63 of Benzing book
% t = parameter taken from pag. 63 of Benzing book

% endplates correction
Lambda_n = Lambda * (1 + 1.9 * h / b);
% 3D correction
Cl_new = Cl / (1 + 2 / Lambda_n);
Cd_new = Cd_0 + (1 / (pi * Lambda_n)) * Cl_new^2;
% rectangular wing correction
Cl_new = Cl_new * (1 + t) ;
Cd_new = Cd_new * (1 + s) ;
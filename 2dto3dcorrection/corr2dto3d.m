function [t,s] = corr2dto3d(lambda)

% CORR2DTO3D calculate t and s correction factors, given lambda
%
% Interpolated functions to calculate t and s.
%
% Example:
%
%   [t,s] = corr2dto3d(lambda)
%
% Reference: fig. 4.02 p. 63. Benzing, Ali-Wings

pt = [-0.000824201967025249,0.0307881216643310,0.0227186958842057];
ps = [-0.000559187516751822,0.0268440121120685,0.0257535758036163];

t = polyval(pt,lambda);
s = polyval(ps,lambda);
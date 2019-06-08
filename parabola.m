function [B] = parabola(alpha_f, Cl_front)

N = length(alpha_f);

XX = [ones(N,1),alpha_f',(alpha_f.^2)'];

YY = Cl_front;

B = XX\YY;

end
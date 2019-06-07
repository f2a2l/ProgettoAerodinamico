function [R] = Rotation (beta)


R = [cos(beta) sin(beta);
    -sin(beta) cos(beta)];

return
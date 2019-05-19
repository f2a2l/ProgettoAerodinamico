function [Cp] = PressureCoeff (v, U)

Cp = 1-(v.^2)./(U^2);

return
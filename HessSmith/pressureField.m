function [cP] = pressureField(uField, vField, U)


[m n] = size(uField);
cP = zeros(m,n);
for i = 1:m
    for j = 1:n
        
        UFieldMag(i,j) = uField(i)^2 + vField(j)^2;
        
        cP(i,j) = 1 - UFieldMag(i,j)/U^2;
    end
end
       
function [stall_or_not,Re] = stall_025(delta_cp, speed, chord)

% STALL is a function for the evaluation of the working condition of an 
% airfoil or a multielement airfoil. Followig the method of Valarezo for
% high lift performance, this function is able to compute whether the
% airfoil is stalling or not.
%
% Output:  
% - stall_or_not       is logical: - true  [1]:  stall
%                                  - false [0]:  NOT stall
%
% Input: 
% - delta_cp           vector containing the pressure coefficients computed
%                      by Hess Smith method
% - speed              speed at which the airfoil is working
% - chord              chord of the airfoil

% Selected Mach number
M = 0.15;

% Compute Reynols number
nu = 1.5e-5;
Re = speed * chord / nu;
Re = round(Re,-5); % Reynols number approximation


% Prandtl - Glauert transformation given the Mach number M
delta_cp = delta_cp / sqrt(1-M^2);

% Compute pressure difference from Hess Smith method
Cp_min = min(delta_cp);
Cp_TE  = delta_cp(end); 
D_cp = abs(Cp_min - Cp_TE);

% Determine maximum pressure difference from Valarezo method for M = 0.15
% using the following empirical data
Reynolds =[1.000  1.100  1.200  1.300  1.400  1.500  1.600  1.700  1.800... 
           1.900  2.000  2.100  2.200  2.300  2.400  2.500  2.600  2.700... 
           2.800  2.900  3.000  3.100  3.200  3.300  3.400  3.500  3.600... 
           3.700  3.800  3.900  4.000  4.100  4.200  4.300  4.400  4.500... 
           4.600  4.700  4.800  4.900  5.000  5.100  5.200  5.300  5.400...
           5.500  5.600  5.700  5.800  5.900  6.000].*1e+06; 

Delta_Cp_max = [7.000 7.124 7.293 7.456 7.668 7.813 7.994 8.133 8.308...
                8.435 8.544 8.695 8.803 8.924 9.045 9.147 9.250 9.377...
                9.468 9.576 9.685 9.776 9.860 9.957 10.05 10.14 10.23...
                10.32 10.40 10.48 10.57 10.68 10.76 10.86 10.95 11.02...
                11.10 11.18 11.26 11.32 11.39 11.45 11.52 11.59 11.66...
                11.71 11.77 11.84 11.89 11.94 12.02];

% Evaluation of maximum delta_Cp from Valarezo method as a function of Re
if Re < Reynolds(1) % if Re is out of the range, the smaller Reynolds numb. 
                    % from Valarezo is considered
    Re = Reynolds(1);
    D_cp_Valarezo = Delta_Cp_max(1);
else
    for i = 1 : length(Reynolds)   
        if Re == Reynolds(i)       
            D_cp_Valarezo = Delta_Cp_max(i);
        end 
    end
end

% Evaluation of stall condition
if D_cp <= D_cp_Valarezo % stall does NOT occur   
    stall_or_not = false;
elseif D_cp > D_cp_Valarezo % stall is occurring
    stall_or_not = true;
end
   

end
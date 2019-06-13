function [stall_or_not,Re] = stall(delta_cp, speed, chord)

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
Reynolds = [1.0000    1.1000    1.2000    1.3000    1.4000    1.5000 ...
            1.6000    1.7000    1.8000    1.9000    2.0000    2.1000 ...
            2.2000    2.3000    2.4000    2.5000    2.6000    2.7000 ...   
            2.8000    2.9000    3.0000    3.1000    3.2000    3.3000 ...
            3.4000    3.5000    3.6000    3.7000    3.8000    3.9000 ...
            4.0000    4.1000    4.2000    4.3000    4.4000    4.5000 ...
            4.6000    4.7000    4.8000    4.9000    5.0000    5.1000 ...   
            5.2000    5.3000    5.4000    5.5000    5.6000    5.7000 ...   
            5.8000    5.9000    6.0000].*1e+6;

Delta_Cp_max = [7.0000    7.9779    8.6156    9.0559    9.4735    9.8347...
               10.2298   10.5572   10.9522   11.1893   11.4376   11.6070...
               11.7988   11.9456   12.0810   12.2729   12.3745   12.5213...
               12.6398   12.7583   12.8881   13.0010   13.1026   13.1816...
               13.2606   13.3340   13.4017   13.4695   13.5428   13.5880...
               13.6331   13.6839   13.7178   13.7686   13.8025   13.8476...
               13.8702   13.8871   13.9041   13.9323   13.9492   13.9548...
               13.9661   13.9774   13.9774   13.9831   13.9831   13.9831...
               13.9774   13.9774   13.9774];

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
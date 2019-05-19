function [t1, t2, t3, t4, t5] = getThickParam(Xt, T, rho_bar, beta_bar)

    % conversion to parameters needed
    rho = rho_bar*(T/Xt)^2;
    t1 = sqrt(2*rho);
    beta = beta_bar * atan(T/(1-Xt));
    
    % generation of matrix
  
    % A = [   Xt,     Xt^2,     Xt^3,     Xt^4;
    %          1,     2*Xt,   3*Xt^2,   4*Xt^3;
    %        0.5,        1,      1.5,        2;
    %          1,        1,        1,        1    ];

    invA =  [               -(2*(2*Xt - 1))/(- Xt^4 + 3*Xt^3 - 3*Xt^2 + Xt),            -1/(Xt^2 - 2*Xt + 1),         -(2*Xt^2)/(Xt^2 - 2*Xt + 1),   -(2*(- Xt^3 + 2*Xt^2))/(Xt^3 - 3*Xt^2 + 3*Xt - 1);
                  -(- 8*Xt^2 + Xt + 1)/(Xt*(- Xt^4 + 3*Xt^3 - 3*Xt^2 + Xt)), (2*Xt + 1)/(Xt^3 - 2*Xt^2 + Xt), (2*(Xt^2 + 2*Xt))/(Xt^2 - 2*Xt + 1),    -(Xt^3 + Xt^2 - 8*Xt)/(Xt^3 - 3*Xt^2 + 3*Xt - 1);
              -(2*(2*Xt^2 + 2*Xt - 1))/(Xt*(- Xt^4 + 3*Xt^3 - 3*Xt^2 + Xt)),  -(Xt + 2)/(Xt^3 - 2*Xt^2 + Xt),   -(2*(2*Xt + 1))/(Xt^2 - 2*Xt + 1), -(2*(- Xt^2 + 2*Xt + 2))/(Xt^3 - 3*Xt^2 + 3*Xt - 1);
                            (3*Xt - 1)/(Xt*(- Xt^4 + 3*Xt^3 - 3*Xt^2 + Xt)),          1/(Xt^3 - 2*Xt^2 + Xt),                 2/(Xt^2 - 2*Xt + 1),                -(Xt - 3)/(Xt^3 - 3*Xt^2 + 3*Xt - 1)];
    
    b = [   T - Xt^0.5*t1;  -0.5*t1*Xt^(-0.5);  -tan(beta/2)-0.25*t1;   -t1     ];
    
    % solution
    t = invA*b;
    
    % unpacking
    t2 = t(1);     
    t3 = t(2); 
    t4 = t(3); 
    t5 = t(4);

end

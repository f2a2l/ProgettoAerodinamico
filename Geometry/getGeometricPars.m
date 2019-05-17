function [gp] = getGeometricPars(param)
%getGeometricPars - returns geometric parameters from profile parameters.
%Syntax:    [gp] = getGeometricPars(param);
%where:     gp = [C, Xc, alpha_te, bc, Xt, T, rho, beta]
%           param = [c1, c2, c3, c4, Xt, T, rho_0_bar, beta_te_bar]

    c1 = param(1);
    c2 = param(2);
    c3 = param(3);
    c4 = param(4);
    Xt = param(5);
    T = param(6);
    rho_bar = param(7);
    beta_bar = param(8);

    keq = @(k) 3*c3*(3*k^2-4*k+1) + 3*c4*(-3*k^2+2*k);
    disp([newline newline 'getGeometricPars/fsolve output'])
    disp('------------------------------')
    Kc = fsolve(keq, 0);
    disp(newline)


    C = 3*c3*Kc*(1-Kc)^2 + 3*c4*(1-Kc)*Kc^2;
    Xc = 3*c3*Kc*(1-Kc)^2 + 3*c4*(1-Kc)*Kc^2;
    alpha_te = atan(c4/(1-c2));
    bc = abs( (6*c3*(3*Kc-2) + 6*c4*(-3*Kc+2)) / (6*c1*(3*Kc-2) + 6*c2*(-3*Kc+2) + 3*Kc^2)^2 );

    rho = rho_bar*(T/Xt)^2;
    beta = beta_bar * atan(T/(1-Xt));

    gp = [C, Xc, alpha_te, bc, Xt, T, rho, beta];

end
function [warnOut, x_transition, Cf] = solverBL(Re, x, y, ue)
%solverBL - inspired by J. Moran's INTGRL; solves integral boundary layer equations
%starting at a stagnation point.
%Uses Thwaites' method for laminar flow, Michel's method to fix transition, Head's method
%for turbulent flow.



    %% external turbulence parameter (as a percentual)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % min: 0.027    --------------->    advised for Re >= 1e6
    % for ncrit = 9: try Tu = 0.0702 (for lower ncrit, raise Tu a bit)
    % for Re = 1e5: try Tu = 0.827 (for best correlation w/ xfoil)
    Tu = 0.075; % [%]



    %% input check
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isscalar(Re)
        error('wrong input; Re must be a scalar number.')
    end
    if ~isvector(x) || ~isvector(y) || ~isvector(ue)
        error('wrong input; either x, y or ue is not a vector.')
    end
    if length(x) ~= length(y)
        error('wrong input; x and y have different lengths.')
    end
    if length(x) ~= length(y)
        error('wrong input; ue does not match mesh dimension.')
    end



    %% define parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx = length(x);


    %% initialisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cf = zeros(1, nx);
    N = Cf;
    warnOut = {};
    xi = getSwiseCoord(x, y); % get streamwise coordinate
    ugrad = gradVel(ue, xi); % get velocity gradient
    ue_fun = @(xxx) interpFun(xi, ue, xxx);
    ugrad_fun = @(xxx) interpFun(xi, ugrad, xxx);



    %% laminar: integration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initial conditions (main variables)
    theta = 0.29234 * sqrt(xi(2)/Re /ue(2));
    % delta = 0.64791 * sqrt(xi(2)/Re /ue(2)); % obtained from value above with Eppler's empirical h(hek)
    h = 0.64791/0.29234;
    eta = 0;

    % initial conditions (auxiliary variables)
    Retheta = Re * ue(2) * theta;

    % trial for initial conditions on eta
    % Ret0 = RethetaCrit(h);
    % deta0 = dn_dret(h);
    % eta = 9 - deta0 * (Ret0 - Retheta);

    % initialisation
    ii = 2; % counter for panels
    x_transition = x(end);
    fopts = optimset('Display','off'); % fsolve options

    % cycle over stations
    while eta < -8.43 - 2.4*log(2.7*tanh(Tu/2.7)/100) 

        % calculate skin friction factor
        Cf(ii) = cflam(Retheta, h);
        N(ii) = eta;

        % go on to next panel, check if such panel exists
        ii = ii + 1;
        if ii > nx
            return
        end

        % integration of bl
        guess_y = [theta; h];
        dxi = xi(ii) - xi(ii-1);

        f = @(n_y) implicitDiffEqn(xi(ii), xi(ii-1), n_y, guess_y, Re, ue_fun, ugrad_fun, 0.9);
        [yy, ~, xflag] = fsolve(f, guess_y, fopts);
        if xflag ~= 1
            warning(['at iteration ' int2str(ii) ', fsolve did not converge (flag ' int2str(xflag) ').'])
        end
        % [yy, ~, errit] = autoLevenMarq(f, guess_y, 10, 1e-6, 1e4, ii);
        % if ~isempty(errit)
        %     disp(['Errors at xi = ' num2str(xi(errit))])
        % end

        % integration of wave amplification
        lambda_o = pgRe(h);
        lambda_n = pgRe(yy(2));
        eta = stepAmplInt(eta, dxi, Retheta, h, theta, lambda_o, Re*yy(1)*ue(ii), yy(2), yy(1), lambda_n, Tu);

        % update variables
        theta = yy(1);
        h = yy(2);
        Retheta = Re * theta * ue(ii);

    end

    % save data
    x_transition = x(ii);



    %% turbulent: data from laminar plate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while ii <= nx

        Rex = Re * xi(ii) * ue(ii);
        Cf(ii) = 0.059 * Rex^(-0.2);

        % update counter
        ii = ii + 1;

    end



    %% postprocessing and returning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cf = Cf .* ue.^2;



end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% calculate streamwise coordinate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xi = getSwiseCoord(x, y)

    nx = length(x);
    
    xi = zeros(1, nx); % streamwise coordinate

    for ii = 2:nx
        dx = x(ii) - x(ii-1);
        dy = y(ii) - y(ii-1);
        xi(ii) = xi(ii-1) + sqrt(dx*dx + dy*dy);
    end
end



%% calculate velocity gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ugrad = gradVel(ue, xi)

    ugrad = zeros(size(ue));

    % first and last step: forward/backward differences (I order)
    ugrad(1) = (ue(2) - ue(1))/(xi(2)-xi(1));
    ugrad(end) = (ue(end) - ue(end-1))/(xi(end) - xi(end-1));

    % all other steps: centered differences (II order)
    for ii = 2:(length(ue)-1)
        ugrad(ii) = (ue(ii+1) - ue(ii-1))/(xi(ii+1) - xi(ii-1));
    end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAMINAR CLOSURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% energy shape parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hs = hek_of_h(h)
    if h < 4.35
        hs = 0.0111*((h-4.35)^2)/(h+1) - 0.0278*((h-4.35)^3)/(h+1) ...
                + 1.528 - 0.0002*((h-4.35)*h)^2;
    else
        hs = 0.015*((h-4.35)^2)/h + 1.528;
    end
end

function dh = dhek_dh(h)
    if h < 4.35
        dh = (111*(2*h - 87/10))/(10000*(h + 1)) - (h^2*(2*h - 87/10))/5000 - (417*(h - 87/20)^2)/(5000*(h + 1)) - (111*(h - 87/20)^2)/(10000*(h + 1)^2) + (139*(h - 87/20)^3)/(5000*(h + 1)^2) - (h*(h - 87/20)^2)/2500;
    else
        dh = (3*(2*h - 87/10))/(200*h) - (3*(h - 87/20)^2)/(200*h^2);
    end
end



%% local skin friction coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cf = cflam(Ret, h)
    if h < 5.5
        lhs = 0.0727*((5.5-h)^3)/(h+1) - 0.07;
    else
        lhs = 0.015*(1 - 1/(h-4.5))^2 - 0.07;
    end
    cf = lhs/Ret;
end



%% local dissipation coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cds = cdiss(Ret, hek, h)
    if h < 4
        D1 = 0.207 + 0.00205 * (4-h)^(5.5);
    else
        D1 = 0.207 - 0.0016 * ((h-4)^2) / (1 + 0.02*(h-4)^2);
    end
    D1 = D1/Ret;
    cds = (hek*D1/2);
end



function dy = lamRHS(xi, y, Re, uefun, ugradfun)

    dy = zeros(size(y)); % preallocation

    % unpack input
    theta = y(1); % theta_(n+1)
    h = y(2); % h_(n+1)

    % get ue, ugrad
    ue = uefun(xi);
    ugrad = ugradfun(xi);

    % auxiliary variables
    Ret = Re * ue * theta; % Re_theta
    hek = hek_of_h(h); % hek(h)

    % momentum thickness derivative from momentum thickness equation
    dy(1) = cflam(Ret,h)/2 - (2+h)*(theta/ue)*ugrad;

    % displacement thickness derivative from energy shape parameter equation
    dy(2) = 2*cdiss(Ret,hek,h) - hek*cflam(Ret,h)/2 - hek*(1-h)*(theta/ue)*ugrad;

end



function f = implicitDiffEqn(xi, xi_o, y, y_o, Re, uefun, ugradfun, biask)

    f = zeros(size(y)); % preallocation

    % unpack input
    dxi = xi - xi_o;
    theta = y(1); % theta_(n+1)
    h = y(2); % h_(n+1)
    theta_o = y_o(1);
    h_o = y_o(2);
    hek_o = hek_of_h(h_o);

    % calculate forward derivative
    d_theta = (theta - theta_o) / dxi;
    d_h = (h - h_o) / dxi;

    % get rhs
    rhs = lamRHS(xi, y, Re, uefun, ugradfun);
    rhs_o = lamRHS(xi_o, y_o, Re, uefun, ugradfun);

    f(1) = d_theta - rhs(1)/2 - rhs_o(1)/2;
    f(2) = (biask*theta*dhek_dh(h) + (1-biask)*theta_o*dhek_dh(h_o)) * d_h ...
            - biask*rhs(2) - (1-biask)*rhs_o(2);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WAVE AMPLIFICATION CLOSURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% derivative of eta with respect to Re theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deriv = dn_dret(h)
    deriv = 0.028*(h-1) - 0.0345*exp(-(3.87/(h-1)-2.52)^2);
end



%% derivative of Retheta with respect to streamwise coordinate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deriv = dret_dx(h)
    deriv = -0.05 + 2.7/(h-1) - 5.5/(h-1)^2 + 3/(h-1)^3;
end



%% critical Re theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret0 = RethetaCrit(h)
    lg10 = (0.267659/(h-1) + 0.394429) * tanh(12.7886/(h-1) - 8.57463) + 3.04212/(h-1) + 0.6660931;
    ret0 = 10^lg10;
end



%% numerical smoothness correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rnorm = RNORM(Ret, h)
    p = 0.08;
    RetCRIT = RethetaCrit(h);
    rnorm = (log10(Ret) - (log10(RetCRIT) - p))/(2*p);
end

function rfac = RFAC(Ret, h)
    rnorm = RNORM(Ret, h);
    if rnorm <=0
        rfac = 0;
    elseif rnorm < 1
        rfac = 3*rnorm^2 - 2*rnorm^3;
    else
        rfac = 1;
    end
end



%% correction for turbulence in freestream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = gturb(Ret, RetONSET)
    A = 0.1;
    B = 0.3;
    r = (Ret/RetONSET - 1)/B + 0.5;
    if r <= 0
        g = 0;
    elseif r < 1
        g = A*(3*r^2 - 2*r^3);
    else
        g = A;
    end
end



%% Retheta onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function re = OnsetRet(lambda, Tu)
    if Tu <= 1.3
        re = (1173.51 - 589.428*Tu + 0.2196/Tu^2)*Flam(lambda, Tu);
    else
        re = 331.5*(Tu - 0.5658)^(-0.671)*Flam(lambda, Tu);
    end
    if re < 20
        re = 20;
    end
end

function F = Flam(lambda, Tu)
    if lambda <= 0
        F = 1 - (-12.986*lambda - 123.66*lambda^2 - 405.689*lambda^3)*exp(-((Tu/1.5)^1.5));
    else
        F = 1 + 0.275*(1-exp(-35*lambda))*exp(-2*Tu);
    end
end

function lambda = pgRe(h)
    lambda = 0.058*(h-4)^2/(h-1) - 0.068;
    if lambda < -0.1
        lambda = -0.1;
        warning('bla') % FIXME: delete me
    elseif lambda > 0.1
        lambda = 0.1;
        warning('bla') % FIXME: delete me
    end
end




%% wave amplification derivative and integration step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rhs = waRHS(Ret, h, theta, lambda, Tu)
    RetONSET = OnsetRet(lambda, Tu);
    rhs = dn_dret(h)*dret_dx(h)*RFAC(Ret, h) + gturb(Ret, RetONSET)/(theta); % if not working, try and comment RFAC
end

function eta = stepAmplInt(eta_o, dxi, Ret_o, h_o, theta_o, lambda_o, Ret, h, theta, lambda, Tu)

    rhs = waRHS(Ret, h, theta, lambda, Tu);
    rhs_o = waRHS(Ret_o, h_o, theta_o, lambda_o, Tu);

    the_real_rhs = rhs/2 + rhs_o/2;

    eta = eta_o + dxi*the_real_rhs;

end
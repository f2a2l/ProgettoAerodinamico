function [warnOut, x_transition, Cf] = solverBL(Re, x, y, ue, varargin)
%solverBL - inspired by J. Moran's INTGRL; solves integral boundary layer equations
%starting at a stagnation point.
%Uses Thwaites' method for laminar flow, Michel's method to fix transition, Head's method
%for turbulent flow.

    h_trans = 1.35;
    L = Re;
    Re_orig = Re;
    Re = 1;
    Re = Re*1e0;

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



    %% load settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % possibly fixed transition
    if isempty(varargin) || varargin{1} == 0
        fixedTrans = false;
    else
        fixedTrans = true;
        xtrans = varargin{1};
    end

    % warning behaviour
    showWarn = true;
    if length(varargin) > 1
        showWarn = varargin{2};
    end


    %% define parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx = length(x);


    %% initialisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cf = zeros(1, nx);
    warnOut = {};
    xi = getSwiseCoord(x, y); % get streamwise coordinate
    ugrad = gradVel(ue, xi); % get velocity gradient



    %% laminar: integration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initial conditions (main variables)
    theta = 0.29004 * sqrt(xi(2)/Re_orig /ue(2)) * L;
    delta = 2.23641 * theta; % derived starting from Eppler's cond. on hek
    eta = 0;

    % initial conditions (auxiliary variables)
    h = 2.23641; % obtained from value above with Eppler's empirical h(hek)
    Retheta = Re * ue(2) * theta;

    % trial for initial conditions on eta
    Ret0 = RethetaCrit(h);
    deta0 = detadre(h);
    eta = 9 - deta0 * (Ret0 - Retheta);
    
    % start stability check
    Rex = Re_orig * xi(2) * ue(2);
    Retmax = 1.174 * (1 + 22400/Rex) * Rex^(0.46);

    % initialisation
    ii = 2; % counter for panels
    x_transition = x(end);

    % cycle over stations
    while eta < 9 && Retheta < Retmax

        % calculate skin friction factor
        Cf(ii) = cflam(Retheta, h, L);

        % go on to next panel, check if such panel exists
        ii = ii + 1;
        if ii > nx
            return
        end

        % integration of bl
        guess_y = [theta; delta];
        dxi = xi(ii) - xi(ii-1);

        % LEVENBERG MARQUARDT
        minthick = min(abs(theta), abs(delta));
        TOL = min(1e-2, minthick) * 1e-4;
        kbias = 1;
        lambda_lq = 10; % 0 for Newton Raphson
        f = @(a) stepLamIntBias(a, Re, theta, delta, dxi, ue(ii-1), ue(ii), ugrad(ii-1), ugrad(ii), kbias, L);
        J = @(a) jacobLamIntBias(a, Re, theta, delta, dxi, ue(ii-1), ue(ii), ugrad(ii-1), ugrad(ii), kbias, L);
        yy = levenMarq(f, J, guess_y, lambda_lq, TOL, 10000, ii);

        % integration of wave amplification
        eta = stepAmplInt(eta, dxi, theta, h, yy(1), yy(2));

        % update variables
        theta = yy(1);
        delta = yy(2);
        h = delta/theta;
        Retheta = Re * theta * ue(ii);

        % check transition
        if fixedTrans
            if x(ii) > xtrans
                break
            end
        else
            Rex = Re_orig * xi(ii) * ue(ii);
            Retmax = 1.174 * (1 + 22400/Rex) * Rex^(0.46);
        end

    end

    % save data
    x_transition = x(ii);
    Cf = Cf/L;



    %% turbulent: data from laminar plate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while ii <= nx

        Rex = Re_orig * xi(ii) * ue(ii);
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
    if h < 4
        hs = 1.515 + 0.076*((4-h)^2)/h;
    else
        hs = 1.515 + 0.040*((h-4)^2)/h;
    end
end

function dh = dhek_dh(h)
    if h < 4
        dh = 0.076*(1-(16/h^2));
    else
        dh = 0.040*(1-(16/h^2));
    end
end



%% local skin friction coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cf = cflam(Ret, h, L)
    if h < 7.4
        lhs = -0.067 + 0.01977*((7.4-h)^2)/(h-1);
    else
        lhs = -0.067 + 0.022*(1 - 1.4/(h-6))^2;
    end
    cf = L * 2 * lhs /Ret;
end

function dcf = dcf_dret(Ret, h, L)
    if h < 7.4
        lhs = -0.067 + 0.01977*((7.4-h)^2)/(h-1);
    else
        lhs = -0.067 + 0.022*(1 - 1.4/(h-6))^2;
    end
    dcf = L * (-2)*lhs/Ret^2;
end

function dcf = dcf_dh(Ret, h, L)
    if h < 7.4
        lhs = 0.01977*((h^2 - 2*h - 39.96)/(h-1)^2);
    else
        lhs = 0.022*(2.8*h - 20.72)/(h-6)^3;
    end
    dcf = L * 2 * lhs /Ret;
end



%% local dissipation coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cds = cdiss(Ret, hek, h, L)
    if h < 4
        lhs = 0.207 + 0.00205 * (4-h)^(5.5);
    else
        lhs = 0.207 - 0.003 * ((h-4)^2) / (1 + 0.02*(h-4)^2);
    end
    cds = L * lhs * hek /2 /Ret;
end

function dcds = dcdiss_dret(Ret, hek, h, L)
    if h < 4
        lhs = 0.207 + 0.00205 * (4-h)^(5.5);
    else
        lhs = 0.207 - 0.003 * ((h-4)^2) / (1 + 0.02*(h-4)^2);
    end
    dcds = L * (-lhs) * hek /2 /Ret^2;
end

function dcds = dcdiss_dhek(Ret, hek, h, L)
    if h < 4
        lhs = 0.207 + 0.00205 * (4-h)^(5.5);
    else
        lhs = 0.207 - 0.003 * ((h-4)^2) / (1 + 0.02*(h-4)^2);
    end
    dcds = L * lhs /2 /Ret;
end

function dcds = dcdiss_dh(Ret, hek, h, L)
    if h < 4
        lhs = 0.00205 * (-5.5)*(4-h)^(4.5);
    else
        lhs = - 0.003 * (5000*h - 20000)/(h^2 - 8*h + 66)^2;
    end
    dcds = L * lhs * hek /2 /Ret;
end



%% biased implicit equation for laminar integration step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = stepLamIntBias(x, Re, theta, delta, dxi, ue, n_ue, ugrad, n_ugrad, kbias, L)

    y = zeros(size(x)); % preallocation

    % LEGEND: plain variables = old variables; n_variables = new variables

    % unpack input
    n_theta = x(1); % theta_(n+1)
    n_delta = x(2); % delta_(n+1)

    % preprocessing: new auxiliary variables
    n_Ret = Re * n_ue * n_theta; % Re_theta
    n_h = n_delta/n_theta; % new h = d/theta
    n_hek = hek_of_h(n_h); % new hek(h)

    % preprocessing: old auxiliary variables
    h = delta/theta; % old h
    hek = hek_of_h(h); % old hek
    Ret = Re * ue * theta;

    % momentum thickness equation: CRANK NICHOLSON
    y(1) = (n_theta - theta)/dxi ...
            + (2+n_h)*(n_theta/n_ue)*n_ugrad/2 - cflam(n_Ret, n_h, L)/4 ...
            + (2+  h)*  (theta/  ue)*  ugrad/2 - cflam(  Ret,   h, L)/4;

    % energy shape parameter equation: BACKWARDS EULER
    y(2) = n_theta*(n_hek-hek)/dxi ...
            + kbias * ((n_hek*(1-n_h))*(n_theta/n_ue)*n_ugrad - 2*cdiss(n_Ret, n_hek, n_h, L) ...
                        + n_hek*cflam(n_Ret, n_h, L)/2) ...
            + (1-kbias) * n_theta/theta*((hek*(1-h))*(theta/ue)*ugrad - 2*cdiss(Ret, hek, h, L) ...
                        + hek*cflam(Ret, h, L)/2);
    
end



%% biased derivatives of equation for laminar integration step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = jacobLamIntBias(x, Re, theta, delta, dxi, ue, n_ue, ugrad, n_ugrad, kbias, L)
    % LEGEND: plain variables = old variables; n_variables = new variables
    
        % preallocation
        J = zeros(2,2);
    
        % unpack input
        n_theta = x(1); % theta_(n+1)
        n_delta = x(2); % delta_(n+1)
    
        % preprocessing: new auxiliary variables
        n_Ret = Re * n_ue * n_theta; % Re_theta
        n_h = n_delta/n_theta; % new h = d/theta
        n_hek = hek_of_h(n_h); % new hek(h)
    
        % preprocessing: old auxiliary variables
        h = delta/theta; % old h
        hek = hek_of_h(h); % old hek
        Ret = Re * ue * theta;
    
        % derivative of momentum thickness equation wrt momentum thickness
        J(1,1) = 1/dxi + n_ugrad/n_ue - dcf_dret(n_Ret,n_h,L)*Re*n_ue/4 + dcf_dh(n_Ret,n_h,L)*n_h/n_theta/4;
    
        % derivative of momentum thickness equation wrt displacement thickness
        J(1,2) = n_ugrad/2/n_ue - dcf_dh(n_Ret,n_h,L)/4/n_theta;
    
        % derivative of energy shape parameter eqn wrt momentum thickness
        J(2,1) = (n_hek - hek)/dxi - n_h/dxi*dhek_dh(n_h) ...
                    + kbias * ( ...
                        + dhek_dh(n_h) * (2*dcdiss_dhek(n_Ret,n_hek,n_h, L)*n_h/n_Ret  ...
                                            - n_ugrad*n_h*(1-n_h)/n_ue - n_h*cflam(n_Ret,n_h,L)/2/n_theta) ...
                        + n_hek * (n_ugrad/n_ue - dcf_dh(n_Ret,n_h,L)*n_h/n_theta/2) ...
                        + Re * n_ue * (n_hek * dcf_dret(n_Ret,n_h,L)/2 - 2*dcdiss_dret(n_Ret,n_hek,n_h, L)) ...
                        + 2 * dcdiss_dh(n_Ret,n_hek,n_h,L) * n_h / n_theta ...
                    ) ...
                    + (1-kbias) /theta*((hek*(1-h))*(theta/ue)*ugrad - 2*cdiss(Ret, hek, h, L) ...
                        + hek*cflam(Ret, h, L)/2);

        % derivative of momentum thickness equation wrt displacement thickness
        J(2,2) = dhek_dh(n_h)/dxi ...
                + kbias * ( ...
                    + dhek_dh(n_h) * ((1-n_h)*n_ugrad/n_ue) - n_hek/n_ue*n_ugrad ...
                    + ((dhek_dh(n_h)*cflam(n_Ret,n_h,L) + n_hek*dcf_dh(n_Ret,n_h,L))/2 ...
                        - 2*(dcdiss_dhek(n_Ret,n_hek,n_h,L)*dhek_dh(n_h) + dcdiss_dh(n_Ret,n_hek,n_h,L)))/n_theta ...
                );

    
    end
    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WAVE AMPLIFICATION CLOSURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% l and m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [l, m] = lm(h)
    l = (6.54*h - 14.07)/h^2;
    m = (0.058*((h-4)^2)/(h-1) - 0.068)/l;
end



%% derivative of eta with respect to Re theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diff = detadre(h)
    diff = 0.01*sqrt(((2.4*h - 3.7 + 2.5*tanh(1.5*h-4.65))^2)+0.25);
end



%% critical Re theta (possilby deprecated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret0 = RethetaCrit(h)
    lhs = (1.415/(h-1)-0.489)*tanh(20/(h-1)-12.9) + 3.295/(h-1) + 0.44;
    ret0 = 10^lhs;
end



%% wave amplification integration step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n_eta = stepAmplInt(eta, dxi, theta, h, n_theta, n_delta)

    n_h = n_delta/n_theta;

    n_d_eta = detadre(n_h);
    [n_l, n_m] = lm(n_h);

    d_eta = detadre(h);
    [l, m] = lm(h);

    n_eta = eta + dxi*( ...
                n_d_eta * (n_m+1)/4 * n_l/n_theta ...
                + d_eta * (  m+1)/4 *   l/  theta ...
            );

end
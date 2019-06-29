function [x_transition, Cf] = solverBL(Re, x, y, ue)
%solverBL - inspired by J. Moran's INTGRL; solves integral boundary layer equations
%starting at a stagnation point.
%Uses Thwaites' method for laminar flow, Michel's method to fix transition, Head's method
%for turbulent flow.



    xi = getSwiseCoord(x, y); % get streamwise coordinate
    ue_fun = @(xixi) interpFun(xi, ue, xixi);





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% USER SETTINGS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %% numerical stabilisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    magicNo = 0.1; % the higher, the more accurate but less stable
    % try 0.94 to get reattachment but not really stable, try 0.1 for best stability
    stab_thresh = 2.5; % value of h after which (h > h_threshold) stabilisation is activated
    dext_fun = @(xxx) 1.72*sqrt(xxx/(Re*ue_fun(xxx))); % Blasius laminar BL
    artReatt = true; % flag for artificial reattachment


    %% numerical integration scheme
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % USING biask <= 1/2 WILL RESULT IN AN UNSTABLE ALGORITHM
    biask = 1; % bias for theta method (0 <= biask <= 1) for energy shape factor eqn.



    %% warning behaviour
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    showWarnings = false; % casts warnings from fsolve
    debugPlot = false; % plots main quantities for debugging





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PROGRAM START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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



    %% pre-allocation and parameters definition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx = length(x); % mesh size
    Cf = zeros(1, nx);
    h_save = Cf;
    theta_save = Cf;
    bret = Cf;
    bretmax = Cf;
    luegrad = gradVel(ue, xi); % get velocity gradient
    luegrad_fun = @(xxx) interpFun(xi, luegrad, xxx);
    CH = channHeight(xi, ue, dext_fun, magicNo); % get fictitious channel height for stabilization



    %% laminar: integration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initial conditions (main variables)
    theta = 0.29234 * sqrt(xi(2)/Re /ue(2));
    h = 0.64791/0.29234;

    % initial conditions (auxiliary variables)
    Retheta = Re * ue(2) * theta;
    Retmax = Retheta + 1;

    % initialisation
    ii = 2; % counter for panels
    x_transition = x(end);
    fopts = optimset('Display','off'); % fsolve options

    % cycle over stations
    while Retheta < Retmax % transition criterion

        % calculate skin friction factor and save previous data
        Cf(ii) = cflam(Retheta, h);
        h_save(ii) = h;
        theta_save(ii) = theta;
        bret(ii) = Retheta;
        bretmax(ii) = Retmax;

        % artificial reattachment
        if artReatt && Cf(ii-1) <0 && Cf(ii) * Cf(ii-1) < 0
            break
        end

        % go on to next panel, check if such panel exists
        ii = ii + 1;
        if ii > nx
            break
        end

        if h < stab_thresh % integration without stabilisation (for better accuracy)

            guess_y = [theta; h];
            f = @(n_y) stepLamInt(xi(ii), xi(ii-1), n_y, guess_y, ue(ii), ue(ii-1), Re, luegrad_fun, biask);
            [yy, ~, xflag] = fsolve(f, guess_y, fopts);
            if showWarnings && xflag ~= 1
                warning(['at iteration ' int2str(ii) ', fsolve did not converge (flag ' int2str(xflag) ').'])
            end

        else % integration with stabilisation

            guess_y = [theta; h; ue(ii-1)];
            f = @(nyy) stepLamInt_wStab(ii, xi, nyy, guess_y, ue, h_save, theta_save, CH, Re);
            [yy, ~, xflag] = fsolve(f, guess_y, fopts);
            if showWarnings && xflag ~= 1
                warning(['at iteration ' int2str(ii) ', fsolve did not converge (flag ' int2str(xflag) ').'])
            end

        end

        % update variables
        if h > stab_thresh
            ue(ii) = yy(3);
        end
        theta = yy(1);
        h = yy(2);
        Retheta = Re * theta * ue(ii);

        % transition criterion
        Rex = Re * xi(ii) * ue(ii);
        Retmax = 1.174 * (1 + 22400/Rex) * Rex^(0.46);

    end

    % save transition coordinate
    if ii <= nx
        x_transition = x(ii);
    end

    % plots for debug
    if debugPlot

        figure
        plot(xi, h_save);
        legend({'h'}) 

        figure
        plot(xi, theta_save);
        legend({'theta'})

        figure
        plot(xi, bret)
        hold on
        plot(xi, bretmax)
        hold off
        legend({'Re theta' 'max Re theta'})
    
    end



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
    ugrad(1) = log(ue(2)/ue(1)) / log(xi(2)/xi(1));
    ugrad(end) = log(ue(end)/ue(end-1)) / log(xi(end)/xi(end-1));

    % all other steps: centered differences (II order)
    for ii = 2:(length(ue)-1)
        ugrad(ii) = log(ue(ii+1)/ue(ii-1)) / log(xi(ii+1)/xi(ii-1));
    end

end



%% calculate fictitious channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ch = channHeight(xi, ue, dext_fun, dotv)
    ch = zeros(size(ue));
    for ii = 1:length(ch)
        d = dext_fun(xi(ii));
        ch(ii) = dotv/ue(ii) + d;
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNSTABILIZED INTEGRATION STEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% right hand side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy = lamRHS(xi, y, Re, ue, luegrad_fun)

    dy = zeros(size(y)); % preallocation

    % unpack input
    theta = y(1); % theta_(n+1)
    h = y(2); % h_(n+1)

    % auxiliary variables
    Ret = Re * ue * theta; % Re_theta
    hek = hek_of_h(h); % hek(h)

    luegrad = luegrad_fun(xi);

    % momentum thickness derivative from momentum thickness equation
    dy(1) = (xi/theta)*(cflam(Ret,h)/2) - (2+h)*luegrad;

    % displacement thickness derivative from energy shape parameter equation
    dy(2) = (xi/theta)*(2*cdiss(Ret,hek,h)/hek - cflam(Ret,h)/2) - (1-h)*luegrad;

end



%% integration step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = stepLamInt(xi, xi_o, y, y_o, ue, ue_o, Re, luegrad_fun, biask)

    f = zeros(size(y)); % preallocation

    % unpack input
    dlxi = log(xi/xi_o);
    theta = y(1); % theta_(n+1)
    h = y(2); % h_(n+1)
    theta_o = y_o(1);
    h_o = y_o(2);

    % auxiliary variables
    hek = hek_of_h(h);
    hek_o = hek_of_h(h_o);

    % calculate forward derivative
    dl_theta = log(theta/theta_o) / dlxi;
    dl_hek = log(hek/hek_o) / dlxi;

    % get rhs
    rhs = lamRHS(xi, y, Re, ue, luegrad_fun);
    rhs_o = lamRHS(xi_o, y_o, Re, ue_o, luegrad_fun);

    f(1) = dl_theta - biask*rhs(1) - (1-biask)*rhs_o(1)/2;
    f(2) = dl_hek ...
           - biask*rhs(2) - (1-biask)*rhs_o(2);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STABILIZED INTEGRATION STEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% right hand side
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy = lamRHS_wStab(xi, y, Re, CH, xi_o, ue_o, h_o, CH_o, theta_o)

    dy = zeros(size(y)); % preallocation

    % unpack input
    theta = y(1); % theta_(n+1)
    h = y(2); % h_(n+1)
    ue = y(3);

    % auxiliary variables
    Ret = Re * ue * theta; % Re_theta
    delta = h * theta;
    hek = hek_of_h(h); % hek(h)
    luegrad = log(ue/ue_o) / log(xi/xi_o); % d(ln(ue))/d(ln(xi))
    lCHgrad = log(CH/CH_o) / log(xi/xi_o); % d(ln(CH))/d(ln(xi))

    % old auxiliary variables
    delta_o = h_o * theta_o;
    ldgrad = log(delta/delta_o) / log(xi/xi_o); % d(ln(delta))/d(ln(xi))
    

    % momentum thickness derivative from momentum thickness equation
    dy(1) = (xi/theta)*(cflam(Ret,h)/2) - (2+h)*luegrad;

    % displacement thickness derivative from energy shape parameter equation
    dy(2) = (xi/theta)*(2*cdiss(Ret,hek,h)/hek - cflam(Ret,h)/2) - (1-h)*luegrad;

    % fictitious ue equation
    dy(3) = (delta*ldgrad - CH*lCHgrad) / (CH-delta);

end



%% integration step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = stepLamInt_wStab(ii, xi_vec, y, y_o, ue_vec, h_vec, theta_vec, CH_vec, Re)

    biask = 1/2;

    f = zeros(size(y)); % preallocation

    % unpack input
    dlxi = log(xi_vec(ii)/xi_vec(ii-1));

    theta = y(1); % theta_(n+1)
    h = y(2); % h_(n+1)
    ue = y(3);
    
    theta_o = y_o(1);
    h_o = y_o(2);
    ue_o = y_o(3);

    % auxiliary variables
    hek = hek_of_h(h);
    hek_o = hek_of_h(h_o);

    % calculate forward derivative
    dl_theta = log(theta/theta_o) / dlxi;
    dl_hek = log(hek/hek_o) / dlxi;
    dl_ue = log(ue/ue_o) / dlxi;

    % get rhs
    rhs = lamRHS_wStab(xi_vec(ii), y, Re, CH_vec(ii), xi_vec(ii-1), ue_vec(ii-1), h_vec(ii-1), CH_vec(ii-1), theta_vec(ii-1));
    rhs_o = lamRHS_wStab(xi_vec(ii-1), y_o, Re, CH_vec(ii-1), xi_vec(ii-2), ue_vec(ii-2), h_vec(ii-2), CH_vec(ii-2), theta_vec(ii-2));

    f(1) = dl_theta - rhs(1)/2 - rhs_o(1)/2;
    f(2) = dl_hek ...
           - biask*rhs(2) - (1-biask)*rhs_o(2);
    f(3) = dl_ue - rhs(3)/2 - rhs_o(3)/2;

end
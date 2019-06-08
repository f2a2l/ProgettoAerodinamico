%function [x_transition] = solverBL(Re, x, y, ue, h_trans, varargin)
function [x_transition, Cf, xi] = solverBL(Re, varargin)
%solverBL - inspired by J. Moran's INTGRL; solves integral boundary layer equations
%starting at a stagnation point.
%Uses Thwaites' method for laminar flow, Michel's method to fix transition, Head's method
%for turbulent flow.

    % TODO: best fit for h_trans; remember it must be between 1.3 and 1.4

    %% input check
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if ~isscalar(Re)
    %     error('wrong input; Re must be a scalar number.')
    % end
    % if ~isvector(x) || ~isvector(y) || ~isvector(ue)
    %     error('wrong input; either x, y or ue is not a vector.')
    % end
    % if length(x) ~= length(y)
    %     error('wrong input; x and y have different lengths.')
    % end
    % if length(x) ~= length(y)
    %     error('wrong input; ue does not match mesh dimension.')
    % end



    %% load settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % possibly fixed transition
    if isempty(varargin)
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
    % nx = length(x); FIXME: debug
    nx = 100;



    %% initialisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = zeros(1, nx);
    Cf = theta;
    warnOut = {};
    % xi = getSwiseCoord(x, y); % get streamwise coordinate
    xi = getXi(nx); % FIXME: debug
    % ugrad = gradVel(ue, xi); % get velocity gradient
    ugrad = getUgrad(xi); %FIXME: debug



    %% laminar: integration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialise
    theta(1) = sqrt(.075/Re/ugrad(1));
    ii = 1; % counter for panels
    Retheta = 0;
    Retheta_max = 1;

    while Retheta < Retheta_max

        % get lambda
        lambda = theta(ii)*theta(ii) * ugrad(ii) * Re;

        % check laminar separation
        if lambda < -0.842
            if showWarn
                warning(['laminar separation at x = ' num2str(x(ii))]);
            end
            warnOut{end+1} = 'laminar separation';
            warnOut{end+1} = x(ii);
            return
        end

        % thwaites' closure
        [l, h] = thwaites(lambda);

        % calculate skin friction factor
        Cf(ii) = 2*l/Re/theta(ii);
        if ii > 1
            Cf(ii) = Cf(ii)/ue(ii);
        end

        % go on to next panel, check if integration is done
        ii = ii + 1;
        if ii > nx
            return
        end

        % integration of bl
        dth2ve6 = .225*(ue(ii)^5 + ue(ii-1)^5) * (xi(ii) - xi(ii-1))/Re;
        theta(ii) = sqrt(((theta(ii-1)*theta(ii-1))*(ue(ii-1)^6) + dth2ve6) / (ue(ii)^6));
        if ii == 2
            theta(2) = theta(1);
        end

        % check transition
        if fixedTrans
            if x(ii) > xtrans
                break
            end
        else
            Rex = Re * xi(ii) * ue(ii);
            Retheta = Re * theta(ii) * ue(ii);
            Retheta_max = 1.174 * (1 + 22400/Rex) * Rex^.46;
        end

    end

    % save data
    x_transition = x(ii);
    %i_trans = ii;



    %% turbulent: integration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialisation
    % h = h_trans;

    disp(['h at transition: ' num2str(h)])
    if h < 1
        warning('you probably need to set h_trans')
    end

    yy(2) = h1_of_h(h);
    yy(1) = theta(ii - 1);

    % keep on integrating
    while ii <= nx

        % get dx and perform integration step
        dx = xi(ii) - xi(ii-1);
        yy = runge2(ii-1, ii, dx, yy, 2, ugrad, Re);

        % unpack integration output
        theta(ii) = yy(1);
        h = h_of_h1(yy(2));

        % get Re and Cf
        Retheta = Re * ue(ii) * theta(ii);
        Cf(ii) = cfturb(Retheta, h);

        % check separation
        if h > 2.4
            if showWarn
                warning(['turbulent separation at x = ' num2str(x(ii))]);
            end
            warnOut{end+1} = 'turbulent separation';
            warnOut{end+1} = x(ii);
            return
        end

    end

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

    for ii = 2:length(x)
        dx = x(ii) - x(ii-1);
        dy = y(ii) - y(ii-1);
        xi(ii) = xi(ii-1) + sqrt(dx*dx + dy*dy);
    end
end



%% differentiate velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  ugrad = gradVel(ue, xi)

    nx = length(ue);

    ugrad = zeros(1, nx); % external velocity gradient

    v1 = ue(3);
    x1 = xi(3);
    v2 = ue(1);
    x2 = xi(1);
    xi(nx+1) = xi(nx-2);
    for ii = 1:nx
        v3 = v1;
        x3 = x1;
        v1 = v2;
        x1 = x2;
        v2 = ue(nx-2);
        if ii < nx
            v2 = ue(ii+1);
        end
        x2 = xi(ii+1);
        fact = (x3 - x1)/(x2 - x1);
        ugrad(ii) = ((v2 - v1) * fact - (v3 - v1)/fact)/(x3-x2);
    end

end



%% laminar closure: thwaites'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l,h] = thwaites(lambda)
    
    if lambda < 0
        l = .22 + 1.402*lambda + .018*lambda/(.107 + lambda);
        h = 2.088 + .0731/(.14 + lambda);
    else
        l = .22 + lambda*(1.57 - 1.8*lambda);
        h = 2.61 - lambda*(3.75 - 5.24*lambda);
    end

end



%% turbulent closure: Head's correlation formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h1 = h1_of_h(h)

    if h > 1.6
        h1 = 3.3 + 1.5501*(h - .6778)^(-3.064);
    else
        h1 = 3.3 + .8234*(h - 1.1)^(-1.287);
    end

end



%% turbulent closure: inverse Head's correlation formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function h = h_of_h1(h1)

    if h1 < 3.3
        h = 3.0;
    elseif h1 < 5.3
        h = .6778 + 1.1536*(h1 - 3.3)^(-.326);
    else
        h = 1.1 + .86*(h1 - 3.3)^(.777);
    end

end



%% turbulent closure: Ludwieg-Tillman skin friction formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function cf = cfturb(rtheta, h)

    cf = .246 * (10 ^(-.678*h)) * rtheta^(-.268);

end



%% turbulent integration: derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function yp = derivs(ii, yt, ugrad, Re)

    h1 = yt(2);

    if h1 <= 3
        error('h1 must be > 3.')
    end

    h = h_of_h1(h1);
    rtheta = Re * ue(ii) * yt(1);

    yp(1) = -(h+2)*yt(1)*ugrad(ii)/ue(ii) + .5*cfturb(rtheta,h);
    yp(2) = -h1*(ugrad(ii)/ue(ii) + yp(1)/yt(1)) + .0306*(h1-3)^(-.6169)/yt(1);

end



%% turbulent integration: Runge-Kutta 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function yy_out = runge2(i0, i1, dx, yy, n, ugrad, Re)

    intvls = i1 - i0;
    
    if intvls < 1
        error('invalid i1 and i0.')
    end

    for ii = 1:intvls
        for jj = 1:n
            yt(jj) = yy(jj);
        end
        yp = derivs(i0+ii-1, yt, ugrad, Re);
        for jj = 1:n 
            yt(jj) = yy(jj) + dx * yp(jj);
            ys(jj) = yy(jj) + .5*dx*yp(jj); 
        end
        yp = derivs(i0+1, yt, ugrad, Re);
        for jj = 1:n
            yy_out(jj) = ys(jj) + .5*dx*yp(jj);
        end
    end

end



function xx = x(ii)
    xx = -cos(pi*(ii-1)/(100-1));
end

function yy = y(ii)
    yy = sin(pi*(ii-1)/100-1) *0.5;
end

function ve = ue(ii)
    ve = (1 + 0.5)*sqrt((1-x(ii)^2)/(1-(1-0.5^2)*x(ii)^2));
end

function xi = getXi(nx)    
    xi = zeros(1, nx); % streamwise coordinate

    for ii = 2:nx
        dx = x(ii) - x(ii-1);
        dy = y(ii) - y(ii-1);
        xi(ii) = xi(ii-1) + sqrt(dx*dx + dy*dy);
    end
end

function ugrad = getUgrad(xi)

    nx = length(xi);
    ugrad = zeros(1, nx); % external velocity gradient

    v1 = ue(3);
    x1 = xi(3);
    v2 = ue(1);
    x2 = xi(1);
    xi(nx+1) = xi(nx-2);
    for ii = 1:nx
        v3 = v1;
        x3 = x1;
        v1 = v2;
        x1 = x2;
        v2 = ue(nx-2);
        if ii < nx
            v2 = ue(ii+1);
        end
        x2 = xi(ii+1);
        fact = (x3 - x1)/(x2 - x1);
        ugrad(ii) = ((v2 - v1) * fact - (v3 - v1)/fact)/(x3-x2);
    end

end
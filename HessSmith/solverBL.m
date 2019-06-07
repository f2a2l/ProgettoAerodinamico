function solverBL(Re, x, y, ue, tau, settings)
%solverBL - inspired by J. Moran's INTGRL; solves integral boundary layer equations
%starting at a stagnation point.
%Uses Thwaites' method for laminar flow, Michel's method to fix transition, Head's method
%for turbulent flow.



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
    if settings{1} == 0
        fixedTrans = false;
    else
        fixedTrans = true;
        xtrans = settings{1};
    end

    % warning behaviour
    showWarn = settings{2};



    %% define parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx = length(x);



    %% initialisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = zeros(1, nx);
    Cf = theta;
    warnOut = {};
    xi = getSwiseCoord(x, y); % get streamwise coordinate
    ugrad = gradVel(ue, xi); % get velocity gradient



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

    
    x_transition = x(ii);


    % TODO: integration of turbulent BL


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
    v2 = ve(1);
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
        ugrad = ((v2 - v1) * fact - (v3 - v1)/fact)/(x3-x2);
    end

end



%% thwaites' laminar flow closure
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
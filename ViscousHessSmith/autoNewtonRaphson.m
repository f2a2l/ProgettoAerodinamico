function  [x,ii] = autoNewtonRaphson(f, x0, TOLL, MAXITER, varargin)

    if ~isempty(varargin)
        dbgFLAG = varargin{1};
    else
        dbgFLAG = false;
    end

    delta = TOLL * 1e-4; % delta_x used to calculate Jacobian at given point

    ii = 1;
    dx = TOLL + 1;

    while dx > TOLL

        J = autoJacob(f, delta, x0);
        dx = J\(-f(x0));
        x = x0 + dx;

        ii = ii + 1;
        if ii > MAXITER
            warning('maximum iterations number reached; algorithm did not converge.')
            return
        end

        dx = norm(dx);
        if dbgFLAG
            gss(ii) = dx;
        end
        x0 = x;

    end

    if dbgFLAG
        figure('Name','autoNewtonRaphson: dx at each iteration')
        plot(1:length(gss), gss)
    end

end



function J = autoJacob(f, delta, x0)

    N = length(f(x0));

    for ii = 1:N
        for jj = 1:N 

            xp = x0;
            xm = x0;
            xp(jj) = xp(jj) + delta;
            xm(jj) = xm(jj) - delta;

            yp = f(xp);
            ym = f(xm);

            J(ii, jj) = (yp(ii) - ym(ii))/2/delta;

        end
    end

    
end
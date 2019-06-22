function  [x,ii] = newtonRaphson(f, df, x0, TOLL, MAXITER)

    ii = 1;
    dx = TOLL + 1;

    while dx > TOLL

        J = df(x0);
        dx = J\(-f(x0));
        x = x0 + dx;

        ii = ii + 1;
        if ii > MAXITER
            warning('maximum iterations number reached; algorithm did not converge.')
            return
        end

        dx = norm(dx);
        x0 = x;

    end

end
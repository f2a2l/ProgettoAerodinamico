function  [x,ii] = newtonRaphson(f, df, x0, TOLL, MAXITER)

    N = length(f(x0));

    ii = 1;
    ff = TOLL + 1;

    lambda = 10;

    while ff > TOLL

        J = df(x0);

        dx = (J'*J + lambda*eye(N)) \ (-J'*f(x0));

        x = x0 + dx;

        ii = ii + 1;
        if ii > MAXITER
            warning('maximum iterations number reached; algorithm did not converge.')
            return
        end

        dx = norm(dx);
        x0 = x;
        ff = norm(f(x0));


    end

    norm(f(x0))
    ii

end
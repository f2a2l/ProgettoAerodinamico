function  [x,ii] = levenMarq(f, df, x0, lambda, TOLL, MAXITER, varargin)

    if ~isempty(varargin)
        outerIter = varargin{1};
    else
        outerIter = [];
    end

    N = length(f(x0));

    ii = 1;
    ff = TOLL + 1;

    while ff > TOLL

        J = df(x0);

        dx = (J'*J + lambda*eye(N)) \ (-J'*f(x0));

        x = x0 + dx;

        ii = ii + 1;
        if ii > MAXITER
            warning(['maximum local iterations number reached at global iteration ' int2str(outerIter) '; algorithm did not converge.'])
            break
        end

        dx = norm(dx);
        x0 = x;
        ff = norm(f(x0));


    end
    
end
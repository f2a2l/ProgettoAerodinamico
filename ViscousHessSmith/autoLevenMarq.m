function  [x,ii] = autoLevenMarq(f, x0, lambda, TOLL, MAXITER, varargin)

    if ~isempty(varargin)
        outerIter = varargin{1};
    else
        outerIter = [];
    end

    N = length(f(x0));

    ii = 1;
    ff = TOLL + 1;

    while ff > TOLL

        J = autoJacob(f, TOLL*1e-4, x0);

        dx = (J'*J + lambda*eye(N)) \ (-J'*f(x0));

        x = x0 + dx;

        ii = ii + 1;
        if ii > MAXITER
            warning(['maximum local iterations number reached at global iteration ' int2str(outerIter) '; algorithm did not converge.'])
            break
        end

        x0 = x;
        ff = norm(f(x0));

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
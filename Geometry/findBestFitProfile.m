function p = findBestFitProfile(x_target, y_target, initial_guess)
    % baseArfl = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
    baseArfl = initial_guess;
    J = @(x) getErr(x, x_target, y_target);
    optz = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e4);
    p = fminunc(J, baseArfl, optz);
    % p = particleswarm(J,8);
    [x,y] = AirfoilShape(p, 80);
    figure('Name', 'findBestFitProfile')
    hold on
    plot(x,y)
    plot(x_target, y_target)
    legend({'IGP', 'target'})
    axis equal
end


function cerr = getErr(p,x_target,y_target)
    [x,y] = AirfoilShape(p, 80);
    cerr = 0;
    for ii = 1:length(x)
        [~, idx] = min(sqrt((x_target - x(ii)).^2 + (y_target - y(ii)).^2));
        cerr = cerr + sqrt( (x(ii) - x_target(idx))^2 + (y(ii) - y_target(idx))^2);
    end
end
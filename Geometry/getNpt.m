function npt = getNpt(orig_npt, relChord)

    % settings
    min_npt = 0.2 * orig_npt;
    fctr = 1.4;

    % get no points
    npt = fctr * relChord * orig_npt;
    npt = floor(npt);
    npt = min([npt orig_npt]);
    npt = max([npt min_npt]);

end
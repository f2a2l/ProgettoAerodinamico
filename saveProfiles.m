if exitflag == 1
    
n_points = 80;

arflPar = [c1_main, c2_main, c3_main, c4_main, x_te_main, T_main, rho_main, beta_te_main;
            c1_flap, c2_flap, c3_flap, c4_flap, x_te_flap, T_flap, rho_flap, beta_te_flap];

profile = solverVHS(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);
profile_DRS = solverVHS_DRS(n_points, arflPar, [AoA_main AoA_flap],[x_flap y_flap],c_flap);

end

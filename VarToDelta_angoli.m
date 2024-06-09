function Delta_max_deg = VarToDelta_angoli(Var_rad2)
    % Var è in rad^2
    Delta_max_deg = (180 * 3 * sqrt(Var_rad2)) / (pi);
    % Delta è in gradi
end
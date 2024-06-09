function Var_rad2 = DeltaToVar_angoli(Delta_max_deg)
    % Delta è in gradi
    Var_rad2 = ( (pi * Delta_max_deg)/(180*3) )^2;
    % Var è in rad^2
end
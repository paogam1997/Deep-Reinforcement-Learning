function quaternion = eulerToQuaternion(DeltaMax_deg)
    % Estrazione dati
    roll_deg = DeltaMax_deg(1);
    pitch_deg = DeltaMax_deg(2);
    yaw_deg = DeltaMax_deg(3);

    % Converti gli angoli di Eulero in radianti
    roll = roll_deg * pi/180;
    pitch = pitch_deg * pi/180;
    yaw = yaw_deg * pi/180;

    % Calcola i valori trigonometrici degli angoli di Eulero
    cy = cos(yaw * 0.5);
    sy = sin(yaw * 0.5);
    cp = cos(pitch * 0.5);
    sp = sin(pitch * 0.5);
    cr = cos(roll * 0.5);
    sr = sin(roll * 0.5);

    % Calcola il quaternione corrispondente
    w = cy * cp * cr + sy * sp * sr;
    x = cy * cp * sr - sy * sp * cr;
    y = sy * cp * sr + cy * sp * cr;
    z = sy * cp * cr - cy * sp * sr;

    quaternion = [w, x, y, z];
end
function euler = quaternionToEuler(w, x, y, z)
    % Calcola gli angoli di Eulero corrispondenti
    roll = atan2(2 * (w * x + y * z), 1 - 2 * (x^2 + y^2));
    pitch = asin(2 * (w * y - z * x));
    yaw = atan2(2 * (w * z + x * y), 1 - 2 * (y^2 + z^2));

    % Converte gli angoli di Eulero in gradi
    roll = rad2deg(roll);
    pitch = rad2deg(pitch);
    yaw = rad2deg(yaw);

    euler = [roll, pitch, yaw];
end
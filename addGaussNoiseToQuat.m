function perturbed_quaternion = addGaussNoiseToQuat(nominal_quaternion, quaternion_var, if_worst_case)
if nargin < 3 % non calcolo per il worst case
    if_worst_case = false;
end
% Check che il quaternione nominale sia normalizzato
if norm(nominal_quaternion) <= 1 + eps || norm(nominal_quaternion) >= 1 - eps
    if if_worst_case == true
        max_quat_variation = 3 * sqrt(quaternion_var);
        perturbed_quaternion_ = nominal_quaternion + max_quat_variation;
        perturbed_quaternion = perturbed_quaternion_ / norm(perturbed_quaternion_);
    else
        % Generiamo un rumore gaussiano per ciascuna componenete del
        % quaternione
        gaussian_noise = sqrt(quaternion_var) .* randn(1, 4);
        % Aggiunge rumore al quaternione nominale
        perturbed_quaternion_ = nominal_quaternion + gaussian_noise;
        % Normalizza il quaternione perturbato
        perturbed_quaternion = perturbed_quaternion_ / norm(perturbed_quaternion_);
    end
else
    disp(" Il quaternione nominale dato non è normalizzato ")
end
end

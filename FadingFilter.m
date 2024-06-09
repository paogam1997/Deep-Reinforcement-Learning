function [x_hat, dx_hat, ddx_hat, dddx_hat] = FadingFilter(beta, previous_x_hat, measured_x, filter_order, varargin)
% Crea un oggetto inputParser
p = inputParser;
% Aggiungi argomenti richiesti
addRequired(p, 'beta');
addRequired(p, 'previous_x_hat');
addRequired(p, 'measured_x');
addRequired(p, 'filter_order');
% Aggiungi argomenti opzionali con valori predefiniti
addOptional(p, 'Ts', 0);
% Analizza gli argomenti di input
parse(p, beta, previous_x_hat, measured_x, filter_order, varargin{:});
% Estrai i valori dagli argomenti di input
Ts = p.Results.Ts;

if filter_order == 1
    % Ridefinizione ingressi per semplicità
    x_hat_p = previous_x_hat(1,:);
    x_star = measured_x(1,:);
    % Definizione del parametro del filtro
    G = 1 - beta;
    % Applicazoine del Fading Filter
    res = x_star - x_hat_p;
    x_hat = x_hat_p + G * res;
    dx_hat = "Non stimabile";
elseif filter_order == 2
    % Ridefinizione ingressi per semplicità
    x_hat_p = previous_x_hat(1,:);
    dx_hat_p = previous_x_hat(2,:);
    x_star = measured_x(1,:);
    % Definizione dei parametri del filtro
    G  = 1 - beta^2;
    H = (1 - beta)^2;
    % Applicazione del Fading Filter del 2^ ordine
    res = x_star - (x_hat_p + dx_hat_p * Ts);
    x_hat = (x_hat_p + dx_hat_p * Ts) + G * res;
    dx_hat = (dx_hat_p) + H/Ts * res;
elseif filter_order == 4
    % Ridefinizione ingressi per semplicità
    x_hat_p = previous_x_hat(1,:);
    dx_hat_p = previous_x_hat(2,:);
    ddx_hat_p = previous_x_hat(3,:);
    dddx_hat_p = previous_x_hat(4,:);
    x_star = measured_x(1,:);
    % Definizione dei parametri del filtro
    G  = beta(1); K = beta(3);
    H = beta(2);  J = beta(4);
    % Applicazione del Fading Filter del 2^ ordine
    res = x_star - (x_hat_p + dx_hat_p*Ts + 0.5*ddx_hat_p*Ts^2 + 1/6 *dddx_hat_p*Ts^3);
    x_hat = (x_hat_p + dx_hat_p*Ts + 0.5*ddx_hat_p*Ts^2 + 1/6 *dddx_hat_p*Ts^3) + G *(res);
    dx_hat = (dx_hat_p + ddx_hat_p*Ts + 0.5*dddx_hat_p*Ts^2) + H/Ts * (res);
    ddx_hat = (ddx_hat_p + dddx_hat_p*Ts) + 2*K/Ts^2 * (res);
    dddx_hat = (dddx_hat_p) + 6*J/Ts^3 * (res);
end
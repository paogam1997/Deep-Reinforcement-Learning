function [x_hat_, P_, e_, svd_L] = KalmanKin(type, dt, Q, R, x_hat_prev, P_prev, y_measured, u_input)
% Applicazione del filtro di Kalman per la stima della cinematica dei
% giunti da una loro misura rumorosa diretta (pulizia del rumore)
% Modello cinematico cambia in base all'input, può essere:
% type = "v_cost" | "a_cost" | "j_cost"
% x_hat_prev è la stima al passo di correzione precedente x(k-1|k-1)
% P_prev è la MSE al passo di correzione precedente P(k-1|k-1)

if strcmp(type, "v_cost")  % Supponiamo modello a velocità costante
    A = [eye(2),     Ts*eye(2);
        zeros(2,2), eye(2)];
    C = eye(4);
    M = eye(4);
    D = [eye(4)];
    %============ Fase di Predizione =============%
    %Calcolo stima predetta al passo successivo applicando il modello della
    %dinamica td (è x_k|k-1)
    x_hat_pred = A * x_hat_prev;
    P_pred = A*P_prev*A' + D*Q*D'; %Matrice MSE della stima al passo k|k-1

    %========== Fase di Correzione ===========%
    %Calcolo la innovazione
    e = y_measured - C * x_hat_pred;
    %    e = [e(1); atan2(sin(e(2)),cos(e(2)))];               %innovazione [rad; rad; rad/s; rad/s]
    S = C*P_pred*C' + M*R*M';                     %matrice di covarianza della innovazione
    L = P_pred*C'*S^(-1);                                      %guadagno di correzione
    x_hat_corr = x_hat_pred + L*e;                                  %calcolo della stima corretta k|k con stima BLUE
    P_corr = (eye(size(A)) - L*C)*P_pred*(eye(size(A)) - L*C)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
    %Salvo i risultati ottenuti ad ogni ciclo per poi graficarli
    x_hat_ = x_hat_corr;
    P_ = P_corr;
    e_ = e;
    svd_L = svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione

elseif strcmp(type, "snap_cost") % Derivata quarta costante
    A = [eye(2), Ts*eye(2), zeros(2,4);
        zeros(2,2), eye(2), Ts*eye(2), zeros(2,2);
        zeros(2,4), eye(2), Ts*eye(2);
        zeros(2, 6), eye(2)];
    C = [eye(4), zeros(4,4)];
    M = eye(4);
    D = [eye(6); zeros(2,6)];

    %============ Fase di Predizione =============%
    %Calcolo stima predetta al passo successivo applicando il modello della
    %dinamica td (è x_k|k-1)
    x_hat_pred = A * x_hat_prev;
    P_pred = A*P_prev*A' + D*Q*D'; %Matrice MSE della stima al passo k|k-1

    %========== Fase di Correzione ===========%
    %Calcolo la innovazione
    e = y_measured - C * x_hat_pred;
    %    e = [e(1); atan2(sin(e(2)),cos(e(2)))];               %innovazione [rad; rad; rad/s; rad/s]
    S = C*P_pred*C' + M*R*M';                     %matrice di covarianza della innovazione
    L = P_pred*C'*S^(-1);                                      %guadagno di correzione
    x_hat_corr = x_hat_pred + L*e;                                  %calcolo della stima corretta k|k con stima BLUE
    P_corr = (eye(size(A)) - L*C)*P_pred*(eye(size(A)) - L*C)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
    %Salvo i risultati ottenuti ad ogni ciclo per poi graficarli
    x_hat_ = x_hat_corr;
    P_ = P_corr;
    e_ = e;
    svd_L = svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione

elseif strcmp(type, "periodic")
    % Modello: dx = [d(somma 3 sin)+w1, d(somma 3 sin)+w2];
    % nel TD: F = I + Ts*d(din)/dx ; D = Ts*d(din)/dw
    F = eye(2);
    D = dt * eye(2);
    H = eye(2);
    M = eye(2);

    %============ Fase di Predizione =============%
    %Calcolo stima predetta al passo successivo applicando il modello della
    %dinamica td (è x_k|k-1)
    x_hat_pred = x_hat_prev + dt * u_input;
    P_pred = F*P_prev*F' + D*Q*D'; %Matrice MSE della stima al passo k|k-1

    %========== Fase di Correzione ===========%
    %Calcolo la innovazione
    e = y_measured - x_hat_pred;
    S = H*P_pred*H' + M*R*M';                     %matrice di covarianza della innovazione
    L = P_pred*H'*S^(-1);                                      %guadagno di correzione
    x_hat_corr = x_hat_pred + L*e;                                  %calcolo della stima corretta k|k con stima BLUE
    P_corr = (eye(size(F)) - L*H)*P_pred*(eye(size(F)) - L*H)' + L*M*R*M'*L';   %Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
    %Salvo i risultati ottenuti ad ogni ciclo per poi graficarli
    x_hat_ = x_hat_corr;
    P_ = P_corr;
    e_ = e;
    svd_L = svd(L);  %salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
end


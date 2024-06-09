%% PARAMETRI CINEMATICA
addpath 'Funzioni Utili'\
clear all
t_max = 3;                       % [s] ampiezza intervallo di simulazione
dt = 0.01;                       % [s] policy step
t = linspace(0,t_max,t_max/dt);  % Genero un campione di tempo per ogni dt da 0 a t_max
%========================= C.I. angoli  =========================%
hip_0 = 120*pi/180;              % [rad] C.I. per hip
knee_0 = 60*pi/180;              % [rad] C.I. per knee

%================ Parametri per imposizione angoli giunti ================%
% Ampiezze scritte come [deg]pi/180 e frequenze come 2pi/[periodo,s]
Hip1 = 15*pi/180;                % [rad] ampiezza oscillazione rollio a BF 
Hip2 = 5*pi/180;                 % [rad] ampiezza oscillazione rollio a MF
Hip3 = 1*pi/180;                 % [rad] ampiezza oscillazione rollio a AF
w_hip1 = 2*pi/2;                 % [rad/s] frequenza oscillazione rollio a BF 
w_hip2 = 2*pi/0.5;               % [rad/s] frequenza oscillazione rollio a MF
w_hip3 = 2*pi/0.2;               % [rad/s] frequenza oscillazione rollio a AF

Knee1 = 20*pi/180;               % [rad] ampiezza oscillazione beccheggio a BF
Knee2 = 15*pi/180;               % [rad] ampiezza oscillazione beccheggio a MF
Knee3 = 3*pi/180;                % [rad] ampiezza oscillazione beccheggio a AF
w_knee1 = 2*pi/1;                % [rad/s] frequenza oscillazione beccheggio a BF
w_knee2 = 2*pi/0.3;              % [rad/s] frequenza oscillazione beccheggio a MF
w_knee3 = 2*pi/0.1;              % [rad/s] frequenza oscillazione beccheggio a AF

%================================= Moto ==================================%
% Moto angoli giunti
hip =   hip_0+   Hip1*sin(w_hip1*t) +     Hip2*sin(w_hip2*t) +     Hip3*sin(w_hip3*t);
knee = knee_0+ Knee1*sin(w_knee1*t) + Knee2*sin(w_knee2*t) + Knee3*sin(w_knee3*t);
% Calcolo della velocità angolare 
dhip = Hip1*w_hip1*cos(w_hip1*t)+Hip2*w_hip2*cos(w_hip2*t)+Hip3*w_hip3*cos(w_hip3*t);
dknee = Knee1*w_knee1*cos(w_knee1*t)+Knee2*w_knee2*cos(w_knee2*t)+Knee3*w_knee3*cos(w_knee3*t);

% Aggiunta del rumore
devstd_hip = 0.25;               % [rad] rumore introdotto dalla policy sul hip joint angle
devstd_knee = 0.25;              % [rad] rumore introdotto dalla policy sul knee joint angle
devstd_dhip = 2.5;               % [rad/s] rumore introdotto dalla policy sulla velocità angolare del hip joint
devstd_dknee = 2.5;              % [rad/s] rumore introdotto dalla policy sulla velocità angolare del knee joint
% Genera il rumore bianco gaussiano a media nulla e variaza specificata per
noise_hip = devstd_hip * randn(size(hip));
noise_knee = devstd_knee * randn(size(knee));
noise_dhip = devstd_dhip * randn(size(dhip));
noise_dknee = devstd_dknee * randn(size(dknee));
% Aggiungi il rumore alla misura
policy_hip = hip + noise_hip;
policy_knee = knee + noise_knee;
policy_dhip = dhip + noise_dhip;
policy_dknee = dknee + noise_dknee;

%% FADING FILTER
%========= Definizione dei parametri di Fading per i vari filtri =========%
% Filtro primo ordine: stimo ogni grandezza singolarmente
beta_hip = 0.9;     % [adim] Fiducia nella stima dell'angolo di hip
beta_knee = 0.9;    % [adim] Fiducia nella dtima dell'angolo di knee
beta_dhip = 0.9;    % [adim] Fiducia nella stima diretta della velocità angolare di hip
beta_dknee = 0.9;   % [adim] Fiducia nella stima diretta della velocità angolare di knee
% Filtro del secondo ordine: Stimo velocità sulla base della lettura dell'angolo
beta_hip2 = 0.85;   % [adim] Fiducia nella stima dell'angolo di hip
beta_knee2 = 0.85;  % [adim] Fiducia nella stima dell'angolo di knee
%==================== Inizializzazione delle variabili ===================%
% Filtro primo ordine
hip_hat = zeros(size(t)); 
hip_hat(1) = hip_0;           % Impongo la prima stima pari alla condizione iniziale che è nota
knee_hat = zeros(size(t)); 
knee_hat(1) = knee_0;         % Impongo la prima stima pari alla condizione iniziale che è nota
dhip_hat = zeros(size(t)); 
dknee_hat = zeros(size(t));
% Filtro del secondo ordine
hip_hat2 = zeros(size(t));     
hip_hat2(1) = hip_0;          % Impongo la prima stima pari alla condizione iniziale che è nota
knee_hat2 = zeros(size(t));    
knee_hat2(1) = knee_0;        % Impongo la prima stima pari alla condizione iniziale che è nota
dhip_hat2 = zeros(size(t));     
dknee_hat2 = zeros(size(t));    
%====================== Implementazione dei filtri =======================%
% Filtro del 1^ ordine
tic                   % Avvio un timer
for i = 2:1:numel(t)  % Itero su tutti gli istanti di tempo
    hip_hat(i) = FadingFilter(beta_hip, hip_hat(i-1), policy_hip(i), 1);         % Stima dell'hip joint angle 
    knee_hat(i) = FadingFilter(beta_knee, knee_hat(i-1), policy_knee(i), 1);     % Stima del knee joint angle
    dhip_hat(i) = FadingFilter(beta_dhip, dhip_hat(i-1), policy_dhip(i), 1);     % Stima dell'hip joint velocity
    dknee_hat(i) = FadingFilter(beta_dknee, dknee_hat(i-1), policy_dknee(i), 1); % Stima della knee joint velocity
end
toc                   % Fermo il timer
% Filtro del 2^ ordine
tic                   % Avvio un timer
for i = 2:1:numel(t)  % Itero su tutti gli istanti di tempo
    [hip_hat2(i), dhip_hat2(i)] = FadingFilter(beta_hip2, [hip_hat2(i-1); dhip_hat2(i-1)], [policy_hip(i); policy_dhip(i)], 2, dt);        % Stima di hip joint angle e velocity
    [knee_hat2(i), dknee_hat2(i)] = FadingFilter(beta_knee2, [knee_hat2(i-1); dknee_hat2(i-1)], [policy_knee(i); policy_dknee(i)], 2, dt); % Stima di knee joint angle e velocity
end
toc                   % Fermo il timer
% PROVA del Filreo del 4^ ordine (non esiste la formula specifica)
% Inizializzazione delle variabili
hip_hat4 = zeros(size(t)); 
dhip_hat4 = zeros(size(t)); 
ddhip_hat4 = zeros(size(t)); 
dddhip_hat4 = zeros(size(t)); 
knee_hat4 = zeros(size(t));  
dknee_hat4 = zeros(size(t));  
ddknee_hat4 = zeros(size(t));  
dddknee_hat4 = zeros(size(t));
% Impongo la prima stima pari alla condizione iniziale nota
hip_hat4(1)  = hip_0; 
knee_hat4(1) = knee_0;
% Definizione dei parametri del filtro
beta4_ = 0; 
G = (1 - beta4_^4);              
H = 1.5*(1-beta4_)^2*(1+beta4_)^2;
K = 0.5*(1-beta4_)^3*(1+beta4_)^1;   
J = 0.5*(1-beta4_)^4;
beta4 = [G, H, K, J];
tic                   % Avvio un timer
for i = 2:1:numel(t)  % Itero su tutti gli istanti di tempo
    [hip_hat4(i), dhip_hat4(i), ddhip_hat4(i), dddhip_hat4(i)] = FadingFilter(beta4, [hip_hat4(i-1); dhip_hat4(i-1); ddhip_hat4(i-1); dddhip_hat4(i-1)], policy_hip(i), 4, 0.01);
    [knee_hat4(i), dknee_hat4(i), ddknee_hat4(i), dddknee_hat4(i)] = FadingFilter(beta4, [knee_hat4(i-1); dknee_hat4(i-1); ddknee_hat4(i-1); dddknee_hat4(i-1)], policy_knee(i), 4, 0.01);
end
toc                   % Fermo il timer

%% KALMAN FILTER
%======================== Definizioni preliminari ========================%
type = "periodic";                        % Tipo di modello da usare: "v_cost" | "snap_cost" | "periodic"
std_dev_dist = 5 * pi/180;                % [rad] disturbo introdotto da inesattezze del modello ed errori di integrazione
Q = (std_dev_dist^2) * eye(2);            % Matrice di covarianza del disturbo nello stato
R = diag([devstd_hip^2, devstd_knee^2]);  % Matrice di covarianza del disturbo di misura
%==================== Inizializzazione delle variabili ===================%
x_hat_res = zeros(2, 1, numel(t));        % Stato stimato
x_hat_KF = [hip_0; knee_0];               % Impongo la prima stima pari alla condizione iniziale nota
P_res = zeros(2, 2, numel(t));            % Matrice di covarianza della stima
P_KF(:,:,1) = Q;                          % C.I. per la matrice di covarianza della stima
e_res = zeros(2, 1, numel(t));            % Residuo di misura
svd_L_res = zeros(2, 1, numel(t));        % Vettore dove raccolgo i valori singolari del guadagno di correzione
detP_res = zeros(1,1,numel(t));           % Vettore dove raccolgo il determinante della matrice di covarianza della stima
% Definizione dei parametri per il modello periodic: vel ang = somma di 3 sinusoidi
% Ampiezze delle sinusoidi
A_hip1 = Hip1*w_hip1;    A_hip2 = Hip2*w_hip2;    A_hip3 = Hip3*w_hip3;
A_knee1 = Knee1*w_knee1;  A_knee2 = Knee2*w_knee2;  A_knee3 = Knee3*w_knee3;
% Frequenze delle sinusoidi
omega_hip1 = w_hip1;    omega_hip2 = w_hip2;    omega_hip3 = w_hip3;
omega_knee1 = w_knee1;  omega_knee2 = w_knee2;  omega_knee3 = w_knee3;
% Sfasamenti delle sinusoidi
psi_hip1 = pi/2;   psi_hip2 = pi/2;   psi_hip3 = pi/2; 
psi_knee1 = pi/2;  psi_knee2 = pi/2;  psi_knee3 = pi/2;

%=================== Applicazione del Filtro di Kalman ===================%
tic   % Avvio del timer
if strcmp(type, "periodic")              % Nel caso di modello periodico facciamo i conti espliciti
    y_meas = [policy_hip; policy_knee];  % Definizione del vettore delle misure
    %u_input = [dhip; dknee];             % Test per debug, uso modello esatto
    % Definizione della forzante in ingresso
    u_input = [A_hip1*sin(omega_hip1*t+psi_hip1)    + A_hip2*sin(omega_hip2*t+psi_hip2)    + A_hip3*sin(omega_hip3*t+psi_hip3);
               A_knee1*sin(omega_knee1*t+psi_knee1) + A_knee2*sin(omega_knee2*t+psi_knee2) + A_knee3*sin(omega_knee3*t+psi_knee3)];
    % Modello: dx = [d(somma 3 sin)+w1, d(somma 3 sin)+w2];
    % nel TD: F = I + Ts*d(din)/dx ; D = Ts*d(din)/dw
    F = eye(2);       % Matrice dinamica nel TD
    D = dt * eye(2);  % Matrice del disturbo nello stato in TD
    % Modello d'osservazione y=h(x,v)=x+v ---> H=dh/dx ; M=dh/dv
    H = eye(2);       
    M = eye(2);
    % Itero per ogni dt
    for  i = 1:1:numel(t)   % Manda avanti indici dei vettori
        if i == 1           % Imposizione della condizione iniziale della stima
            x_hat_KF(:,:,i) = [hip_0; knee_0];
            P_KF(:,:,i) = Q;
            e = 0;
            L = zeros(2,2);
        else
            %==================== Fase di Predizione =====================%
            %Calcolo stima predetta al passo successivo applicando il modello della
            %dinamica td (è x_k|k-1)
            Bu = u_input(:,i);                % Scritto così per comodità nel debug
            x_hat_KF = x_hat_KF + dt * Bu ;   % Applicazione del modello della dinamica TD
            P_KF = F*P_KF*(F') + D*Q*(D');    % Matrice MSE della stima al passo k|k-1
            %==================== Fase di Correzione =====================%
            y_star = y_meas(:,i);             % Scritto così per comodità nel debug
            e =  y_star - x_hat_KF;           % Calcolo la innovazione
            S = H*P_KF*(H') + M*R*(M');       % Matrice di covarianza della innovazione
            L = P_KF*(H')*(S^(-1));           % Guadagno di correzione
            x_hat_KF = x_hat_KF + L*e;        % Calcolo della stima corretta k|k con stima BLUE
            P_KF = (eye(size(F)) - L*H)*P_KF*(eye(size(F)) - L*H)' + L*M*R*(M')*(L');   % Matrice MSE della stima al passo k|k sfruttando la forma di Joseph
        end
        % Salvo i risultati ottenuti ad ogni ciclo per poi graficarli
        x_hat_res(:,:,i) = x_hat_KF;
        P_res(:,:,i) = P_KF;
        detP_res(:,:,i) = det(P_KF);
        e_res(:,:,i) = e;
        svd_L_res(:,:,i) = svd(L);  % salvo i valori singolari del guadagno di correzione per eventuale verifica di effettivo lavoro della fase di correzione
    end
else % Negli altri casi si può usare la funzione appositamente creata
    [x_hat_KF(:,:,i), P_KF(:,:,i), e_KF(:,:,i), svd_L_KF(:,:,i)] = KalmanKin(type, dt, Q, R, x_hat_KF(:,:,i-1), P_KF(:,:,i-1), y_meas(:,i), u_input(:,i));
end
toc     % Fermo il timer
% Estrazione dei segnali stimati
hip_hat_KF = squeeze(x_hat_res(1,1,:))'; 
knee_hat_KF = squeeze(x_hat_res(2,1,:))';

%% PLOT
fig1 = figure(1);
sgtitle("Hip & Knee joint Angles & Velocities")
subplot(2,2,1)
    plot(t, hip, "Color",'r',"LineStyle", "--", "LineWidth", 1);hold on;
    plot(t, hip_hat, "Color", "#0072BD", "LineStyle", "-"); 
    plot(t, hip_hat2, "Color", "#77AC30", "LineStyle", "-");
    plot(t, hip_hat_KF, "Color", "#7E2F8E" )
    hold off; grid on; box on; title("Hip Angle"); xlabel("time [s]"); ylabel("Joint angle [rad]"); legend("Ideale", "FF 1^{st} order", "FF 2^{nd} order", "KF ")
subplot(2,2,3)   
    plot(t, knee, "Color","r","LineStyle","--","LineWidth",1);hold on;  
    plot(t, knee_hat, "Color", "#0072BD", "LineStyle", "-"); 
    plot(t, knee_hat2, "Color", "#77AC30", "LineStyle", "-");
    plot(t, knee_hat_KF, "Color", "#7E2F8E" )

    hold off; grid on; box on; title("Knee Angle"); xlabel("time [s]"); ylabel("Joint angle [rad]"); legend("Ideale", "FF 1^{st} order", "FF 2^{nd} order", "KF ")
subplot(2,2,2)
    plot(t, dhip, "Color",'r',"LineStyle","--","LineWidth",1); hold on;
    plot(t, dhip_hat, "Color", "#0072BD", "LineStyle", "-"); 
    plot(t, dhip_hat2, "Color", "#77AC30", "LineStyle", "-");
    hold off; grid on; box on; title("Hip Velocity"); xlabel("time [s]"); ylabel("Joint velocity [rad/s]"); legend("Ideale", "FF 1^{st} order", "FF 2^{nd} order")
subplot(2,2,4)
    plot(t, dknee, "Color",'r',"LineStyle","--","LineWidth",1); hold on;
    plot(t, dknee_hat, "Color", "#0072BD", "LineStyle", "-"); 
    plot(t, dknee_hat2, "Color", "#77AC30", "LineStyle", "-");
    hold off; grid on; box on; title("Knee Velocity"); xlabel("time [s]"); ylabel("Joint velocity [rad/s]"); legend("Ideale", "FF 1^{st} order", "FF 2^{nd} order")

fig2 = figure(2);
sgtitle("Errori di stima")
subplot(2, 2, 1)
    plot(t, (hip_hat - hip) * 180/pi, "Color","#0072BD"); hold on;
    plot(t, (hip_hat2 - hip) * 180/pi, "Color","#77AC30");
    plot(t, (hip_hat_KF - hip)* 180/pi, "Color", "#7E2F8E" );
    hold off; grid on; box on; title("Hip Angle"); xlabel("time [s]"); ylabel("Joint angle error [deg]"); legend("FF 1^{st} order", "FF 2^{nd} order", "KF ")
subplot(2, 2, 3)
    plot(t, (knee_hat - knee) * 180/pi, "Color","#0072BD"); hold on;
    plot(t, (knee_hat2 - knee) * 180/pi, "Color","#77AC30");
    plot(t, (knee_hat_KF - knee)* 180/pi, "Color", "#7E2F8E" );
    hold off; grid on; box on; title("Knee Angle"); xlabel("time [s]"); ylabel("Joint angle error [deg]"); legend("FF 1^{st} order", "FF 2^{nd} order", "KF ")
subplot(2, 2, 2)
    plot(t, (dhip_hat - dhip) * 180/pi, "Color","#0072BD"); hold on;
    plot(t, (dhip_hat2 - dhip) * 180/pi, "Color","#77AC30");
    hold off; grid on; box on; title("Hip Velocity"); xlabel("time [s]"); ylabel("Joint velocity error [deg/s]"); legend("FF 1^{st} order", "FF 2^{nd} order")
subplot(2, 2, 4)
    plot(t, (dknee_hat - dknee) * 180/pi, "Color","#0072BD"); hold on;
    plot(t, (dknee_hat2 - dknee) * 180/pi, "Color","#77AC30");
    hold off; grid on; box on; title("Knee Velocity"); xlabel("time [s]"); ylabel("Joint velocity error [deg/s]"); legend("FF 1^{st} order", "FF 2^{nd} order")

fig3 = figure(3);
sgtitle("Medie & Deviazioni Standard [deg]-[deg/s]")
subplot(2,2,1)
bar([mean(hip), mean(policy_hip), mean(hip_hat), mean(hip_hat2), mean(hip_hat_KF); ...
     mean(knee), mean(policy_knee), mean(knee_hat), mean(knee_hat2), mean(knee_hat_KF); ...
     mean(dhip), mean(policy_dhip), mean(dhip_hat), mean(dhip_hat2), 0; ...
     mean(dknee), mean(policy_dknee), mean(dknee_hat), mean(dknee_hat2), 0;].*180/pi);
title("Signal Mean");xticklabels({'Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'}) ;legend({'Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord', "KF "},"Location","northeast"); grid on; box on;
subplot(2,2,2)
bar([std(hip), std(policy_hip), std(hip_hat), std(hip_hat2), std(hip_hat_KF);
     std(knee) std(policy_knee) std(knee_hat) std(knee_hat2), std(knee_hat_KF); ...
     std(dhip), std(policy_dhip), std(dhip_hat), std(dhip_hat2), 0; 
     std(dknee) std(policy_dknee) std(dknee_hat) std(dknee_hat2), 0].*180/pi)
title("Signal Standard Deviation");xticklabels({'Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'}) ;legend({'Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord', "KF "},"Location","northwest"); grid on; box on;
subplot(2,2,3)
bar([0, mean(policy_hip - hip), mean(hip_hat - hip), mean(hip_hat2 - hip), mean(hip_hat_KF - hip); ...
     0, mean(policy_knee - knee), mean(knee_hat - knee), mean(knee_hat2 - knee), mean(knee_hat_KF - knee); ...
     0, mean(policy_dhip - dhip), mean(dhip_hat - dhip), mean(dhip_hat2 - dhip), 0; ...
     0, mean(policy_dknee - dknee), mean(dknee_hat - dknee), mean(dknee_hat2 - dknee), 0].*180/pi)
title("Error Mean");xticklabels({'Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'}) ;legend({'Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord', "KF"},"Location","northwest"); grid on; box on;
subplot(2,2,4)
bar([0, std(policy_hip - hip), std(hip_hat - hip), std(hip_hat2 - hip), std(hip_hat_KF - hip); ...
     0 std(policy_knee - knee) std(knee_hat - knee) std(knee_hat2 - knee), std(knee_hat_KF - knee); ...
     0, std(policy_dhip - dhip), std(dhip_hat - dhip), std(dhip_hat2 - dhip), 0; ...
     0 std(policy_dknee - dknee) std(dknee_hat - dknee) std(dknee_hat2 - dknee), 0].*180/pi)
title("Error Standard Deviation");xticklabels({'Hip Angle', 'Knee Angle', 'Hip Vel', 'Knee Vel'}) ;legend({'Ideale', 'Policy', 'FF 1^{st} ord','FF 2^{nd} ord', "KF "},"Location","northwest"); grid on; box on;

%===== Grafico innovazione per controllare sia un rumore bianco =====%
j = 1;
ylabel_ = {'e(1) [rad]'; 'e(2) [rad]'; 'e(3) [rad/s]'; 'e(4) [rad/s]' };
title_ = {'innovazione hip angle'; 'innovazione knee angle'; 'innovazione hip vel'; 'innovazione knee vel'};
for k = 4:(4+size(e_res,1)-1)
    figure(k)
    subplot(2,1,1)
    plot(t,squeeze(e_res(j,1,:)));xlabel('tempo [s]');ylabel(ylabel_{j})
    title(title_{j});
    subplot(2,1,2)
    autocorr(squeeze(e_res(j,1,:)));
    j = j+1;
end

figure(4+size(e_res,1))
subplot(1,2,1); 
plot(t,squeeze(svd_L_res(1,:,:)))  %valore singolare massimo della matrice di guadagno di correzione
title("max svd L"); grid on; box on;
subplot(1,2,2); 
plot(t, squeeze(detP_res))
title("det P"); grid on; box on;

 %% SIMULAZIONE
 prompt="Avviare simulazione doppio pendolo? Y/N [Y]:";
txt=input(prompt,"s");
% Oltre a 'Y' anche premere solo invio è visto come un consenso
if isempty(txt)
    txt='Y';
end
if txt=='Y'
    disp('Ok simulo')

 O = [0;0];
 lhip = 0.175;
 lknee = 0.19;
figure(); hold on; xlim([-0.5,0.5]);xlabel('orizzontale [m]'); ylim([-0.5,0.1]);ylabel('verticale [m]'); grid on;box on; axis equal; title('blu=vero  rosso=policy  verde = FF1 rosa = FF2  nero=FK'); plot(-0.5,-0.5);plot(-0.5,0.1);plot(0.5,-0.5);plot(0.5,0.1);
or = [0; 0];                                         %Definizione della origine           
plot_origine = plot(or(1), or(2), 'ob');             %Plot della origine
 for i=1:1:numel(t)
   [Hip, Knee] = Kin(lhip, lknee, hip(i), knee(i));
   %grafico i due link ed i due giunti
    plot_upperlink = plot([or(1),Hip(1)],[or(2),Hip(2)],'b',"LineWidth",1);
    plot_hipjoint = plot(Hip(1), Hip(2), 'ob',"LineWidth",1);
    plot_lowerlink = plot([Hip(1), Knee(1)], [Hip(2), Knee(2)], 'b',"LineWidth",1);
    plot_feet = plot(Knee(1), Knee(2), '*b',"LineWidth",1);

    %Inserisco anche la traiettoria data da policy
    p = [policy_hip(i), policy_knee(i)];                  %variabili di giunto stimate
    [o1, o2] = Kin(lhip, lknee,p(1), p(2));      %calcolo posizione stimata dei giunti
    %sovrappongo i due link e giunti
    plot_link1_p = plot([or(1), o1(1)],[or(2), o1(2)], 'r--',"LineWidth",1);
    plot_giunto2_p = plot(o1(1), o1(2), 'hr',"LineWidth",1);
    plot_link2_p = plot([o1(1), o2(1)], [o1(2), o2(2)], 'r--',"LineWidth",1);
    plot_ee_p = plot(o2(1), o2(2), 'hr',"LineWidth",1);

    %Inserisco anche la traiettoria data da FF 1^ordine
    ff1 = [hip_hat(i), knee_hat(i)];                  %variabili di giunto stimate
    [o1_ff1, o2_ff1] = Kin(lhip, lknee,ff1(1), ff1(2));      %calcolo posizione stimata dei giunti
    %sovrappongo i due link e giunti
    plot_link1_ff1 = plot([or(1), o1_ff1(1)],[or(2), o1_ff1(2)], 'g--',"LineWidth",1);
    plot_giunto2_ff1 = plot(o1_ff1(1), o1_ff1(2), 'hg',"LineWidth",1);
    plot_link2_ff1 = plot([o1_ff1(1), o2_ff1(1)], [o1_ff1(2), o2_ff1(2)], 'g--',"LineWidth",1);
    plot_ee_ff1 = plot(o2_ff1(1), o2_ff1(2), 'hg',"LineWidth",1);

    %Inserisco anche la traiettoria data da FF 2^ordine
    ff2 = [hip_hat2(i), knee_hat2(i)];                  %variabili di giunto stimate
    [o1_ff2, o2_ff2] = Kin(lhip, lknee,ff2(1), ff2(2));      %calcolo posizione stimata dei giunti
    %sovrappongo i due link e giunti
    plot_link1_ff2 = plot([or(1), o1_ff2(1)],[or(2), o1_ff2(2)], 'm--',"LineWidth",1);
    plot_giunto2_ff2 = plot(o1_ff2(1), o1_ff2(2), 'hm',"LineWidth",1);
    plot_link2_ff2 = plot([o1_ff2(1), o2_ff2(1)], [o1_ff2(2), o2_ff2(2)], 'm--',"LineWidth",1);
    plot_ee_ff2 = plot(o2_ff2(1), o2_ff2(2), 'hm',"LineWidth",1);

    %Inserisco anche la traiettoria data da KF
    kf = [hip_hat_KF(i), knee_hat_KF(i)];                  %variabili di giunto stimate
    [o1_kf, o2_kf] = Kin(lhip, lknee,kf(1), kf(2));      %calcolo posizione stimata dei giunti
    %sovrappongo i due link e giunti
    plot_link1_kf = plot([or(1), o1_kf(1)],[or(2), o1_kf(2)], 'k--',"LineWidth",1);
    plot_giunto2_kf = plot(o1_kf(1), o1_kf(2), 'hk',"LineWidth",1);
    plot_link2_kf = plot([o1_kf(1), o2_kf(1)], [o1_kf(2), o2_kf(2)], 'k--',"LineWidth",1);
    plot_ee_kf = plot(o2_kf(1), o2_kf(2), 'hk',"LineWidth",1);
    timer=text(-0.4,0,"Timer: "+num2str(t(i),2));  %aggiungo il timer

    pause (0.001);                                      %dopo 0.001 secondi cancella i grafici attivi e prosrgui sostituendoli con quelli della iterazione successiva
    %cancella i grafici seguenti
    delete(plot_upperlink);
    delete(plot_hipjoint);
    delete(plot_lowerlink);
    delete(plot_feet);
    delete(plot_link1_p);
    delete(plot_giunto2_p);
    delete(plot_link2_p);
    delete(plot_ee_p);
    delete(plot_link1_ff1);
    delete(plot_giunto2_ff1);
    delete(plot_link2_ff1);
    delete(plot_ee_ff1);
    delete(plot_link1_ff2);
    delete(plot_giunto2_ff2);
    delete(plot_link2_ff2);
    delete(plot_ee_ff2);
    delete(plot_link1_kf);
    delete(plot_giunto2_kf);
    delete(plot_link2_kf);
    delete(plot_ee_kf);
    delete(timer);
 end
hold off
else
     disp('Ok niente simulazione') %messaggio che compare se non simuli
end
function [OHip, OKnee] = Kin(lhip, lknee, hip, knee)
 OHip = lhip * [cos(-hip); sin(-hip)];
 OKnee = OHip + lknee * [cos(-(hip-knee)); sin(-(hip-knee))];
end





% Parametri di comando

vxcomm = rand(1);                                % [m/s] velocità longitudinale desiderata
yawcomm = (2 * rand(1) - 1) * pi;                % [rad] yaw desiderato in [-pi; pi]
dt = 0.005;                                      % [s] passo d'integrazione
yawratecomm = (2 * rand(1) - 1) * 1;       % [rad/s] yaw rate desiderato in [-1; 1]

disp("v_x desiderata"); disp(vxcomm)             % Manda a schermo la v_x desiderata
disp("yaw desiderato"); disp(yawcomm)            % Manda a schermo lo yaw desiderato
disp("yaw_rate desiderato"); disp(yawratecomm)   % Manda a schermo lo yaw_rate desiderato

%% Scale factor della reward

% ==================== 1. Linear Velocity X Tracking ==================== %
linearVelocityXYRewardScale = 1.5;                               % [adim] max reward value
varlinvel =  0.02;                                               % [m^2/s^2] gaissian variance

% ==================== 2. Linear Velocity Z Penalty ===================== %
linearVelocityZRewardScale = -4.0;                               % [adim] quadratic penality coeff
vzRewardScalep = 0.05;                                           % [adim] local gaussian bonus max
fBc_vz = 35;                                                     % [s^2/m^2] penality-bonus blending coeff
varvz = -1/(fBc_vz*linearVelocityZRewardScale);                  % [m^2/s^2] gaussian variance

% =================== 3. Angular Velocity XY Penalty ==================== %
angularVelocityXYRewardScale = -0.05;                            % [adim] quadratic penality coeff
AngularVelocityRewardScalep = 0.05;                              % [adim] local gaussian bonus max
fBc_rollrate = 850;                                              % [s^2/rad^2]  penality-bonus blending coeff
varrollrate = -1/(fBc_rollrate*angularVelocityXYRewardScale);    % [rad^2/s^2] gaussian variance
fBc_pitchrate = fBc_rollrate;                                    % [s^2/rad^2]  penality-bonus blending coeff
varpitchrate = -1/(fBc_rollrate*angularVelocityXYRewardScale);   % [rad^2/s^2] gaussian variance

% =================== 4. Angular Velocity Z Tracking ==================== %
angularVelocityZPenalityScale = -0.75;
angularVelocityZRewardScale = 0.75;                              % [adim] max reward value
fBc_angvel = 150;                                                   % [rad^2/s^2] gaussian variance
varangvel = -1/(fBc_angvel*angularVelocityZPenalityScale);

% ======================= 5. Orientation Penalty ======================== %
orientationRewardScale = -0.75;                                  % [adim] quadratic penality coeff
orientRewScalep = 0.025;                                         % [adim] local gaussian bonus max  
fBc_orient = 50;                                                 % [1/rad^2]  penality-bonus blending coeff
varorient = -1/(fBc_orient*orientationRewardScale);              % [rad^2] gaussian variance

% ========================== 6. Torque Penalty ========================== %
torqueRewardScale = -0.;                                         % [adim] quadratic penality coeff

% ======================= 7. Torque Rate Penalty ======================== %
torque_rateRewardScale = -0.;                                    % [adim] quadratic penality coeff
torque_rateRewardScale_p = 0.;                                   % [adim] local gaussian bonus max 
fBc_torq_rate = 500;                                             % [s^2/Nm^2]  penality-bonus blending coeff
vartorqrate = -1/(fBc_torq_rate*torque_rateRewardScale);         % [Nm^2/s^2] gaussian variance

% ======================= 8. Base Height Penalty ======================== %
baseHeightRewardScale = -10.;                                     % [adim] quadratic penality coeff
baseHeightRewardScale1 = 0.05;                                     % [adim] local gaussian bonus max
fBc_base = 50;                                                   % [1/m^2]  penality-bonus blending coeff
varbase = -1/(fBc_base*baseHeightRewardScale);                   % [m^2] gaussian variance

% ========================= 9. Falcata Penalty ========================== %
falcataRewardScale = -10;                                        % [adim] quadratic penality coeff
falcataRewScalep = 0.1;                                          % [adim] local gaussian bonus max
fBc_falcata = 1;                                                 % [1/m^2]  penality-bonus blending coeff
varfalcata = -1/(fBc_falcata*falcataRewardScale);                % [m^2] gaussian variance

% ======================= 10. Action Rate Penalty ======================= %
actionRateRewardScale = -0.00075;                                                     % [adim] quadratic penality coeff
actRateRewScalep = 0.05;                                                              % [adim] local gaussian bonus max
actRateRewScaleGN = -0.053;                                                           % [adim] local gaussian penality minimum 
fBc_kneeactionRate = 100;                                                             % [1/rad^2]  penality-bonus blending coeff
fBc_scale = 5.2;                                                                      % [adim] scale coeff of local penality-local bonus blending coeff
varkneeActRate = -1/(fBc_kneeactionRate*actionRateRewardScale);                       % [rad^2] gaussian variance
fBc_hipactionRate = fBc_kneeactionRate;                                               % [1/rad^2]  penality-bonus blending coeff
varhipActRate = -1/(fBc_hipactionRate*actionRateRewardScale);                         % [rad^2] local bonus gaussian variance
fBc_actRateGN = -fBc_scale*(actRateRewScaleGN/actRateRewScalep)*fBc_kneeactionRate;   % [1/rad^2]  local penality-local bonus blending coeff
varActRateGN = -1/(fBc_actRateGN*actionRateRewardScale);                              % [rad^2] local penality gaussian variance

% ======================= 11. Joint Accel Penalty ======================= %
jointAccRewardScale = -0.0015;                                   % [adim] quadratic penality coeff  
jointAccRewScalep = 0.05;                                        % [adim] local gaussian bonus max
fBc_kneeAccel = 255;                                             % [s^2/rad^2]  penality-bonus blending coeff
varkneeAccel = -1/(fBc_kneeAccel*jointAccRewardScale);           % [rad^2/s^2] gaussian variance
fBc_hipAccel = fBc_kneeAccel;                                    % [s^2/rad^2]  penality-bonus blending coeff
varhipAccel = -1/(fBc_hipAccel*jointAccRewardScale);             % [rad^2/s^2] gaussian variance

% ======================= 12. Action Accel Penalty ====================== %
actAccRewardScale = -0.07;                                       % [adim] quadratic penality coeff 
actAccRewScalep = 0.1;                                           % [adim] local gaussian bonus max
fBc_kneeActAccel = 70;                                           % [s^2/rad^2]  penality-bonus blending coeff
varkneeActAccel = -1/(fBc_kneeAccel*actAccRewardScale);          % [rad^2/s^2] gaussian variance
fBc_hipActAccel = fBc_kneeAccel;                                 % [s^2/rad^2]  penality-bonus blending coeff
varhipActAccel = -1/(fBc_hipActAccel*actAccRewardScale);         % [rad^2/s^2] gaussian variance

% ======================= 13. Joint Jerk Penalty ========================= %
jointJerkRewardScale = -0.0001;                                  % [adim] quadratic penality coeff 
jointJerkRewardScalep = 0.0375;                                  % [adim] local gaussian bonus max 
fBc_kneej = 20;                                                  % [s^4/rad^2]  penality-bonus blending coeff
varjerkknee = -1/(fBc_kneej*jointJerkRewardScale);               % [rad^2/s^4] gaussian variance
fBc_hipj = fBc_kneej;                                            % [s^4/rad^2]  penality-bonus blending coeff
varjerkhip = -1/(fBc_hipj*jointJerkRewardScale);                 % [rad^2/s^4] gaussian variance

% ======================= 14. Action Jerk Penalty ======================= %
actJerkRewardScale = -0.0005;                                    % [adim] quadratic penality coeff  
actJerkRewardScalep = 0.075;                                     % [adim] local gaussian bonus max
fBc_actkneej = 10;                                               % [s^4/rad^2]  penality-bonus blending coeff
varactjerkknee = -1/(fBc_actkneej*actJerkRewardScale);           % [rad^2/s^4] gaussian variance
fBc_acthipj = fBc_actkneej;                                      % [s^4/rad^2]  penality-bonus blending coeff
varactjerkhip = -1/(fBc_acthipj*actJerkRewardScale);             % [rad^2/s^4] gaussian variance

% ===================== 15. Action Max Jerk Penalty ===================== %    
jerk_maxRewardScale = -0.012;                                    % [adim] quadratic penality coeff   
jerk_maxRewardScale_p = 0;                                       % [adim] local gaussian bonus max
fBc_jerk_max = 500;                                              % [s^4/rad^2]  penality-bonus blending coeff
varjerk_max = -1/(fBc_jerk_max*jerk_maxRewardScale);             % [rad^2/s^4] gaussian variance

% ===================== 16. Joint Max Jerk Penalty ====================== %    
jjerk_maxRewardScale = -0.0018;                                  % [adim] quadratic penality coeff   
jjerk_maxRewardScale_p = 0;                                      % [adim] local gaussian bonus max
fBc_jjerk_max = 500;                                             % [s^4/rad^2]  penality-bonus blending coeff
varjjerk_max = -1/(fBc_jjerk_max*jjerk_maxRewardScale);          % [rad^2/s^4] gaussian variance

% ======================= 17. Feet Pos Z Penalty ======================== %
feetZRewardScale = -40;                                          % [adim] quadratic penality coeff
feetZRewScalep = 0.1;                                            % [adim] local gaussian bonus max
fBc_feetZ = 60;                                                  % [1/m^2]  penality-bonus blending coeff
varfeetZ = -1/(fBc_feetZ*feetZRewardScale);                      % [m^2] gaussian variance

% ======================= 18. Fallen Over Penalty ======================= %
fallenOverRewardScale = -1000.0;                                 % [adim] quadratic penality coeff

% ===================== 19. Action Max Accel Penalty ==================== %    
accel_maxRewardScale = -0.7;                                     % [adim] quadratic penality coeff   
accel_maxRewardScale_p = 0;                                      % [adim] local gaussian bonus max
fBc_accel_max = 500;                                             % [s^2/rad^2]  penality-bonus blending coeff
varaccel_max = -1/(fBc_accel_max*accel_maxRewardScale);          % [rad^2/s^2] gaussian variance

% ===================== 20. Joint Max Accel Penalty ===================== %    
jaccel_maxRewardScale = -0.015;                                  % [adim] quadratic penality coeff   
jaccel_maxRewardScale_p = 0;                                     % [adim] local gaussian bonus max
fBc_jaccel_max = 500;                                            % [s^2/rad^2]  penality-bonus blending coeff
varjaccel_max = -1/(fBc_jaccel_max*jaccel_maxRewardScale);       % [rad^2/s^2] gaussian variance

% ===================== 21. Feet Air Time Penalty ===================== %    
airtime_RewardScale = -0.7;                                       % [adim] quadratic penality coeff   
airtime_RewardScale_p = 0.1;                                     % [adim] local gaussian bonus max
fBc_airtime = 70;                                                 % [1/s^2]  penality-bonus blending coeff
varairtime = -1/(fBc_airtime*airtime_RewardScale);               % [s^2] gaussian variance

%% Parametri di sistema
 
vxlim = 1;                             % [m/s] massima velocità longitudinale
yawmax = pi/6;                         % [rad] massimo valore di yaw accettato
yawratelim = 1;                        % [rad/s] limite per lo yawrate
tau_rate_max = 5;                      % [Nm/s] rate di coppia massima generabile dai motori
initbaseheight = 0.365;                % [m] altezza iniziale base
rollmax = pi/6;                        % [rad] massimo rollio ammissibile
pitchmax = pi/6;                       % [rad] massimo beccheggio ammissibile
kneeaccelmax = pi/9/dt;                % [rad/s] massima velocità angolare relativa per knee tra 2 istanti successivi
hipaccelmax = pi/18/dt;                % [rad/s] massima velocità angolare relativa per hip tra 2 istanti successivi
actionkneemax = 5*pi/3;                % [rad] massima differenza tra action relative al knee tra 2 istanti successivi
actionhipmax = 5*pi/3;                 % [rad] massima differenza tra action relative al hip tra 2 istanti successivi
numrobots = 4096;                      % [adim] numero robot simulati in contemporanea
desbaseheight = 0.316;                 % [m] altezza da terra desiderata per la base
hipjerkmax = (0.1*hipaccelmax)/dt;     % [rad/s^2] massimo jerk sull'hip visto come differenza di accelerazioni angolari
kneejerkmax = (0.1*kneeaccelmax)/dt;   % [rad/s^2] massimo jerk sul knee visto come differenza di accelerazioni angolari
falcata_des = 0.18;                    % [m] ampiezza media della falcata per ogni gamba
feetZ_des = 0.12;                      % [m] altezza z rispetto alla base (positiva se sotto la base) desiderata per il piede a fine falcata 
airtime_des = 0.7;                     % [s] tempo medio desiderato per la falcata

%% Variabili di stato / osservazioni

vx = linspace(-vxlim, vxlim, 500);                                 % [m/s] velocità attuale longitudinale robot
yawrate = linspace(-1,1);                   % [rad/s] yaw rate attuale
% heading = yawrate * dt;                                          % [rad] direzione attuale
% omegaz_comm = 0.5*wrapToPi(yawcomm - heading);                   % yaw rate desiderato
%yawratecomm = max(min(omegaz_comm, yawratelim), -yawratelim);     % yaw rate desiderato limitato tra -1 e 1
tauknee_rate = linspace(-tau_rate_max, tau_rate_max);              % [Nm] coppia sul ginnocchio attuale
tauhip_rate = linspace(-tau_rate_max, tau_rate_max);               % [Nm] coppia sul hip attuale
vz = linspace(-initbaseheight/0.8, initbaseheight/0.8, 10000);     % [m/s] velocità di abbassamento della base
rollrate = linspace(-rollmax/dt, rollmax/dt);                      % [rad/s] roll rate attuale
%rollrate = linspace(-3, 3);
pitchrate = linspace(-pitchmax/dt, pitchmax/dt);                   % [rad/s] pitch rate attuale
%pitchrate = rollrate;
%kneeaccel = linspace(-kneeaccelmax, kneeaccelmax);                % [rad/s] accel angolare knee vista come differenza di velocità angolari!
kneeaccel = linspace(-10, 10);hipaccel = kneeaccel;
% hipaccel = linspace(-hipaccelmax, hipaccelmax);                  % [rad/s] accel angolare hip vista come differenza di velocità angolari!
deltactionknee = linspace(-actionkneemax, actionkneemax,200);      % [rad] velocità angolare knee vista come differenza di action
%deltactionknee = linspace(-0.5,0.5); deltactionhip = deltactionknee;
deltactionhip = linspace(-actionhipmax, actionhipmax,200);         % [rad] velocità angolare hip vista come differenza di action
fallen = linspace(0, numrobots);                                   % [adim] numero di robot caduti
baseheight = linspace(-initbaseheight, initbaseheight);            % [m] altezza attuale della base 
%hipjerk = linspace(-hipjerkmax, hipjerkmax);                      % [rad/s^2] jerk del giunto di hip come diferenza di accelerazioni angolari
hipjerk = linspace(-25, 25); kneejerk = hipjerk;
%kneejerk = linspace(-kneejerkmax, kneejerkmax);                   % [rad/s^2] jerk del giunto di knee come differenza di accelerazioni angolari
prgrav_x = linspace(-1, 1);                                        % [adim] componente x della gravità proiettata nel frame body
prgrav_y = linspace(-1, 1);                                        % [adim] componente y della gravità proiettata nel frame body
falcata = linspace(-0.365, 0.365);                                 % [m] falcata al massimo va pari alla lunghezza della gamba stesa = 0.365
feet_z = linspace(-0.365, 0);                                      % [m] il piede non si alza sopra la base ed al massimo tocca terra -> valore minimo = - altezza base
airtime = linspace(0, 2);                                          % [s] durata di timpo in cui il piede sta in aria 

%% Funzioni di reward (singolo robot)

% ==================== 1. Linear Velocity X Tracking ==================== %
% Ricompensa più la velocità longitudinale della base è vicina alla
% desiderata. Calcolata come solo una gaussiana positiva
rewlinvelxy = linearVelocityXYRewardScale * exp(-(vxcomm - vx).^2 /(2*varlinvel));

% ==================== 2. Linear Velocity Z Penalty ===================== %
% Penalizza moto verticale della base --> problema moto viene penalizzato
% anche se stesse cercando di riportarsi ad una configurazione alta --> va
% in contrasto con il termine di base height
rewlinvelz = linearVelocityZRewardScale * vz.^2 + vzRewardScalep*exp(-vz.^2 /(2*varvz));
% idea di modifica:se la base è abbastanza alta penalizzagli il moto
for k = 1:length(baseheight)
    if baseheight(k) >= 0.8*desbaseheight && baseheight(k) <= 1.2*desbaseheight
        rewlinvelzmod = linearVelocityZRewardScale * vz.^2 + vzRewardScalep*exp(-vz.^2 /(2*varvz));
    else % se la base sta bassa allora non pesare il movimento così riesce a salire
        rewlinvelzmod = 0;
    end
end
clear k;

% =================== 3. Angular Velocity XY Penalty ==================== %
% Penalizza moti di rollio e beccheggio
rewangvelxy = zeros();
for i = 1:length(rollrate)
    for j = 1: length(pitchrate)
        rewangvelxy(i,j) = angularVelocityXYRewardScale * (rollrate(i)^2 + pitchrate(j)^2) + AngularVelocityRewardScalep*exp(-rollrate(i)^2 /(2*varrollrate)) *exp(-pitchrate(j)^2 /(2*varpitchrate));
    end
end
clear i j ;

% =================== 4. Angular Velocity Z Tracking ==================== %
% Ricompensa più la velocità angolare della base è vicina a quella
% desiderata
rewangvelz = angularVelocityZPenalityScale *(yawratecomm - yawrate).^2 + angularVelocityZRewardScale * exp(-(yawratecomm - yawrate).^2 /(2*varangvel));

% ======================= 5. Orientation Penalty ======================== %
reworient = zeros();
for i = 1: length(prgrav_x)
    for j = 1: length(prgrav_y)
        reworient(i,j) = orientationRewardScale * (prgrav_x(i)^2 + prgrav_y(j)^2) + orientRewScalep * exp( -(prgrav_x(i)^2 + prgrav_y(j)^2)/(2*varorient) );
    end
end
clear i j;

% ====================== 6-7. Torque Rate Penalty ======================= %
% Penalizza l'uso di coppie alte su knee e hip
rewtorque_rate = zeros();
for i = 1:length(tauknee_rate)
    for j = 1: length(tauhip_rate)
        rewtorque_rate(i,j) = torque_rateRewardScale_p* exp(-0.5*max(abs(tauknee_rate(i)),abs(tauhip_rate(j))).^2/vartorqrate) + torque_rateRewardScale * max(abs(tauknee_rate(i)) , abs(tauhip_rate(j))) ; % per singola gamba
    end
end
clear i j ;

% ======================= 8. Base Height Penalty ======================== %
% Penalizza il robot più la sua altezza è diversa dalla desiderata
rewbaseheight = baseHeightRewardScale * (baseheight - desbaseheight).^2;
% variante bastone-carota
rewbaseheight1 = baseHeightRewardScale * (baseheight - desbaseheight).^2 + baseHeightRewardScale1 * exp(-(baseheight - desbaseheight).^2 /(2*varbase));

% ========================= 9. Falcata Penalty ========================== %
% Reward per la falcata
errfalcata = (falcata - falcata_des).^2;
rewfalcata = errfalcata * falcataRewardScale + exp(-errfalcata/(2*varfalcata)) * falcataRewScalep;

% ======================= 10. Action Rate Penalty ======================= %
% Penalizza la velocità angolare dei giunti
rewactionrate = zeros();
for i = 1: length(deltactionknee)
    for j = 1: length(deltactionhip)
        rewactionrate(i,j) = actionRateRewardScale * (deltactionknee(i)^2 + deltactionhip(j)^2) + actRateRewScalep*exp(-deltactionknee(i)^2/(2*varkneeActRate)) *exp(-deltactionhip(j)^2/(2*varhipActRate)) + actRateRewScaleGN*exp(-deltactionknee(i)^2/(2*varActRateGN)) *exp(-deltactionhip(j)^2/(2*varActRateGN));
    end
end
rewactionrate_mono = actionRateRewardScale * (deltactionhip.^2) + actRateRewScalep*exp(-deltactionhip.^2/(2*varhipActRate)) + actRateRewScaleGN*exp(-deltactionhip.^2/(2*varActRateGN));
clear i j ;

% ======================= 11. Joint Accel Penalty ======================= %
% Penalizza la accelerazione angolare dei giunti
rewjointaccel = zeros();
for i = 1: length(kneeaccel)
    for j = 1: length(hipaccel)
        rewjointaccel(i,j) = jointAccRewardScale * (kneeaccel(i)^2 + hipaccel(j)^2) + jointAccRewScalep*exp(-kneeaccel(i)^2/(2*varkneeAccel)) *exp(-hipaccel(j)^2/(2*varhipAccel));
    end
end
clear i j ;

% ====================== 12. Action Accel Penalty ======================= %
% Penalizza la accelerazione angolare del moto desiderato dei giunti
rewactaccel = zeros();
for i = 1: length(kneeaccel)
    for j = 1: length(hipaccel)
        rewactaccel(i,j) = actAccRewardScale * (kneeaccel(i)^2 + hipaccel(j)^2) + actAccRewScalep*exp(-kneeaccel(i)^2/(2*varkneeActAccel)) *exp(-hipaccel(j)^2/(2*varhipActAccel));
    end
end
clear i j ;

% ======================= 13. Joint Jerk Penalty ======================== %
% Penalizza il jerk angolare dei giunti
rewjointjerk = zeros();
for i = 1: length(kneejerk)
    for j = 1: length(hipjerk)
        rewjointjerk(i,j) = jointJerkRewardScale * (kneejerk(i)^2 + hipjerk(j)^2) + jointJerkRewardScalep*exp(-kneejerk(i)^2 /(2*varjerkknee)) *exp(-hipjerk(j)^2 /(2*varjerkhip)) ;
    end
end
clear i j ;

% ====================== 14. Action Jerk Penalty ======================== %
% Penalizza il jerk angolare del moto desiderato dei giunti
rewactjerk = zeros();
for i = 1: length(kneejerk)
    for j = 1: length(hipjerk)
        rewactjerk(i,j) = actJerkRewardScale * (kneejerk(i)^2 + hipjerk(j)^2) + actJerkRewardScalep*exp(-kneejerk(i)^2 /(2*varactjerkknee)) *exp(-hipjerk(j)^2 /(2*varactjerkhip)) ;
    end
end
clear i j ;

% ===================== 15. Action Max Jerk Penalty ===================== %
% Penalizza l'uso di picchi di jerk
rewjerk_max = zeros();
for i = 1:length(kneejerk)
    for j = 1: length(hipjerk)
        rewjerk_max(i,j) = jerk_maxRewardScale_p* exp(-0.5*max(abs(kneejerk(i)),abs(hipjerk(j))).^2/varjerk_max) + jerk_maxRewardScale * max(abs(kneejerk(i)),abs(hipjerk(j))) ; % per singola gamba
    end
end
clear i j ;

% ===================== 16. Joint Max Jerk Penalty ====================== %
% Penalizza l'uso di picchi di jerk nel moto attuale
rewjjerk_max = zeros();
for i = 1:length(kneejerk)
    for j = 1: length(hipjerk)
        rewjjerk_max(i,j) = jjerk_maxRewardScale_p* exp(-0.5*max(abs(kneejerk(i)),abs(hipjerk(j))).^2/varjjerk_max) + jjerk_maxRewardScale * max(abs(kneejerk(i)),abs(hipjerk(j))) ; % per singola gamba
    end
end
clear i j ;

% ======================= 17. Feet Pos Z Penalty ======================== %
% Reward per la feet_pos_z
errfeetZ = (feet_z + desbaseheight - feetZ_des).^2;
rewfeetZ = errfeetZ * feetZRewardScale + exp(-errfeetZ/(2*varfeetZ)) * feetZRewScalep;

% ======================= 18. Fallen Over Penalty ======================= %
% Penalizza la caduta del robot
rewfallenover = fallenOverRewardScale * fallen;

% ===================== 19. Action Max Accel Penalty ===================== %
% Penalizza l'uso di picchi di jerk
rewaccel_max = zeros();
for i = 1:length(kneeaccel)
    for j = 1: length(hipaccel)
        rewaccel_max(i,j) = accel_maxRewardScale_p* exp(-0.5*max(abs(kneeaccel(i)),abs(hipaccel(j))).^2/varaccel_max) + accel_maxRewardScale * max(abs(kneeaccel(i)),abs(hipaccel(j))) ; % per singola gamba
    end
end
clear i j ;

% ===================== 20. Joint Max Accel Penalty ====================== %
% Penalizza l'uso di picchi di jerk nel moto attuale
rewjaccel_max = zeros();
for i = 1:length(kneeaccel)
    for j = 1: length(hipaccel)
        rewjaccel_max(i,j) = jaccel_maxRewardScale_p* exp(-0.5*max(abs(kneeaccel(i)),abs(hipaccel(j))).^2/varjaccel_max) + jaccel_maxRewardScale * max(abs(kneeaccel(i)),abs(hipaccel(j))) ; % per singola gamba
    end
end
clear i j ;

% ===================== 21. Feet Air Time Penalty ====================== %
% Penalizza passi troppo veloci
rewairtime = airtime_RewardScale_p* exp(-0.5*(airtime-airtime_des).^2/varairtime) + airtime_RewardScale * (airtime-airtime_des).^2 ; % per singola gamba

% ========================= Utility ===================================== %
[x, y] = meshgrid(-10:1:10, -10:1:10);   % Crea una griglia di punti nel piano x-y
z = zeros(size(x));  % Poiché il piano è parallelo al piano x-y a z=0, tutti i valori di z sono 0

%% Plot

%======================== 1. Termini monodimensionali ====================%
figure(1)
subplot(2, 3, 1); plot(vx, rewlinvelxy); hold on; xline(vxcomm, '--r');
    grid on; box on; title("reward velocità longitudinale"); xlabel("vx [m/s]");xlim([-vxlim, vxlim]); ylabel("rew [adim]"); ylim([0, 1.2*linearVelocityXYRewardScale]); legend('reward', 'vx_{command}'); hold off;

subplot(2,3,2); plot(yawrate, rewangvelz); hold on; xline(yawratecomm, '--r'); 
    grid on; box on; title("reward yaw rate"); xlabel("yawrate [rad/s]"); xlim([-yawratelim, yawratelim]); ylabel("rew [adim]"); ylim([0, 1.1*angularVelocityZRewardScale]); legend('reward', 'yaw rate_{command}'); hold off;

subplot(2,3,3); plot(vz, rewlinvelz); hold on; xline(0, '--r'); %plot(vz, rewlinvelzmod, 'k.--')
    grid on; box on; title("reward velocità verticale"); xlabel("vz [m/s]");xlim([min(vz), max(vz)]); ylabel("rew [adim]"); legend('reward', 'vz_{command}'); hold off;

subplot(2,3,4); plot(fallen, rewfallenover);
    grid on; box on; title("reward fallen over"); xlabel("fallen robots [adim]"); xlim([0, numrobots]); ylabel("rew [adim]"); 

subplot(2,3,5); plot(baseheight, rewbaseheight); hold on; xline(desbaseheight, '--r'); plot(baseheight, rewbaseheight1, 'k.--')
    grid on; box on; title("reward altezza base"); xlabel("z [m]");xlim([0, initbaseheight]); ylabel("rew [adim]"); legend('reward', 'z_{desiderata}', 'mod rew'); hold off;

subplot(2,3,6); plot(falcata, rewfalcata); hold on; xline(falcata_des, '--r');
    grid on; box on; title("reward falcata"); xlabel("falcata [m]"); xlim([min(falcata), max(falcata)]); ylabel("rew [adim]"); legend('reward', 'falc. des.'); hold off
% figure(19)
% plot(feet_z+desbaseheight, rewfeetZ); hold on; xline(feetZ_des, '--r');
%     grid on; box on; title("reward feetZ"); xlabel("altezza piede da terra [m]"); xlim([0, max(feet_z+desbaseheight)]); ylabel("rew [adim]"); legend('reward', 'alt. des.'); hold off
%fig21 = figure(21);
%plot(airtime, rewairtime); hold on; xline(airtime_des, '--r');
%    grid on; box on; title("reward air time"); xlabel("feet air time [s]"), xlim([min(airtime), max(airtime)]); ylabel("rew [adim]"); legend('reward', 'air time des'); hold off

%========================= 2. Reward Torque Rate =========================%
% fig2 = figure(2); 
% sgtitle("Reward Torque Rate")
% subplot(2,2,1); surf(tauknee_rate, tauhip_rate, rewtorque_rate); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
%     box on; grid on; xlabel('\Delta\tau_{i} [Nm]'); xlim([-tau_rate_max, tau_rate_max]);
%     set(gca, 'XDir','reverse'); ylabel('\Delta\tau_{i+1} [Nm]'); ylim([-tau_rate_max, tau_rate_max]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% 
% subplot(2,2,2); pcolor(tauknee_rate, tauhip_rate, rewtorque_rate);
%     box on; grid on; xlabel('\Delta\tau_{i} [Nm]'); xlim([-tau_rate_max, tau_rate_max]); set(gca, 'XDir','reverse'); ylabel('\Delta\tau_{i+1} [Nm]'); ylim([-tau_rate_max, tau_rate_max]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% 
% subplot(2,2,3); plot(tauknee_rate, rewtorque_rate, 'b'); hold on; xline(0, '--r'); yline(0, '--g');
%     box on; grid on; xlabel('\Delta\tau_{i-knee} [Nm]'); xlim([-tau_rate_max, tau_rate_max]); ylabel('rew [adim]'); hold off
% 
% subplot(2,2,4); plot(tauhip_rate, rewtorque_rate,'b'); hold on; xline(0, '--r'); yline(0, '--g');
%     box on; grid on; xlabel('\Delta\tau_{i-hip} [Nm]'); xlim([-tau_rate_max, tau_rate_max]); ylabel('rew [adim]'); hold off
% 
% fontsize(fig2, scale=1.2)   

% %===================== 3. Reward Angular Velocity xy =====================%
% fig3 = figure(3);
% sgtitle("Reward Angular Velocity xy")
% subplot(2,2,1); surf(rollrate, pitchrate, rewangvelxy); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
%     box on; grid on; xlabel('roll rate [rad/s]'); xlim([min(rollrate), max(rollrate)]); set(gca, 'XDir','reverse'); ylabel('pitch rate [rad/s]'); ylim([min(pitchrate), max(pitchrate)]); set(gca,'YDir','reverse'); zlabel('rew [adim]'); hold off
% 
% subplot(2,2,2); pcolor(rollrate, pitchrate, rewangvelxy);
%     box on; grid on; xlabel('roll rate [rad/s]'); xlim([min(rollrate) , max(rollrate)]); set(gca, 'XDir','reverse'); ylabel('pitch rate [rad/s]'); ylim([min(pitchrate), max(pitchrate)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% 
% subplot(2,2,3); plot(rollrate, rewangvelxy, 'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('roll rate [rad/s]'); xlim([min(rollrate) , max(rollrate)]); ylabel('rew [adim]'); hold off
% 
% subplot(2,2,4); plot(pitchrate, rewangvelxy,'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('pitch rate [rad/s]'); xlim([min(pitchrate), max(pitchrate)]); ylabel('rew [adim]'); hold off
% 
% fontsize(fig3, scale=1.2)   
%     
% %================== 4. Reward Joint Angular Acceleration =================%
% fig4 = figure(4);    
% sgtitle("Reward Action Angular Acceleration")
% subplot(2,2,1); surf(kneeaccel, hipaccel, rewactaccel); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
%     box on; grid on; xlabel('\Delta knee_{action rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action rate } [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% 
% subplot(2,2,2); pcolor(kneeaccel, hipaccel, rewactaccel);
%    box on; grid on; xlabel('\Delta knee_{action rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action rate} [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% 
% subplot(2,2,3); plot(kneeaccel, rewactaccel, 'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta knee_{action rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); ylabel('rew [adim]'); hold off
% 
% subplot(2,2,4); plot(hipaccel, rewactaccel,'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta hip_{action rate} [rad/s]'); xlim([min(hipaccel), max(hipaccel)]); ylabel('rew [adim]'); hold off
% 
% fontsize(fig4, scale=1.2)   
% 
% %================== 5. Reward Action Angular Acceleration =================%
% fig5 = figure(5);    
% sgtitle("Reward Joint Angular Acceleration")
% subplot(2,2,1); surf(kneeaccel, hipaccel, rewjointaccel); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
%     box on; grid on; xlabel('\Delta knee_{joint vel} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint vel} [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% 
% subplot(2,2,2); pcolor(kneeaccel, hipaccel, rewjointaccel);
%    box on; grid on; xlabel('\Delta knee_{joint vel} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint vel} [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% 
% subplot(2,2,3); plot(kneeaccel, rewjointaccel, 'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta knee_{joint vel} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); ylabel('rew [adim]'); hold off
% 
% subplot(2,2,4); plot(hipaccel, rewjointaccel,'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta hip_{joint vel} [rad/s]'); xlim([min(hipaccel), max(hipaccel)]); ylabel('rew [adim]'); hold off
% 
% fontsize(fig5, scale=1.2)   
% 
% % %=================== 6. Reward Joint Angular Velocity ====================%
% % fig6 = figure(6);   
% % sgtitle("Reward Joint Angular Velocity")
% % subplot(2,2,1); surf(deltactionknee, deltactionhip, rewactionrate);  hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
% %     box on; grid on; xlabel('\Delta knee_{joint ang} [rad/s]'); xlim([min(deltactionknee) , max(deltactionknee)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint ang} [rad/s]'); ylim([min(deltactionhip), max(deltactionhip)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% % 
% % subplot(2,2,2); pcolor(deltactionknee, deltactionhip, rewactionrate);
% %     box on; grid on; xlabel('\Delta knee_{joint ang} [rad/s]'); xlim([min(deltactionknee) , max(deltactionknee)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint ang} [rad/s]'); ylim([min(deltactionhip), max(deltactionhip)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% % subplot(2,2,3); plot(deltactionknee, rewactionrate, 'b'); hold on; xline(0, '--r'); yline(0, '--g');
% %     box on; grid on; xlabel('\Delta knee_{joint ang} [rad/s]'); xlim([min(deltactionknee) , max(deltactionknee)]); ylabel('rew [adim]'); hold off
% % 
% % subplot(2,2,4); plot(deltactionhip, rewactionrate,'b'); hold on; yline(0, '--g'); xline(0, 'r--');
% %     box on; grid on; xlabel('\Delta hip_{joint ang} [rad/s]'); xlim([min(deltactionhip), max(deltactionhip)]); ylabel('rew [adim]'); hold off
% % 
% % fontsize(fig6, scale=1.2)   
% 
% % %===================== 7. Reward Joint Angular Jerk ======================%
% % fig7 = figure(7);   
% % sgtitle("Reward Joint Angular Jerk")
% % subplot(2,2,1); surf(kneejerk, hipjerk, rewjointjerk); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
% %     box on; grid on; xlabel('\Delta knee_{joint accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint accel} [rad]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% % 
% % subplot(2,2,2); pcolor(kneejerk, hipjerk, rewjointjerk);
% %     box on; grid on; xlabel('\Delta knee_{joint accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint accel} [rad]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% % 
% % subplot(2,2,3); plot(kneejerk, rewjointjerk, 'b'); hold on; yline(0, '--g'); xline(0, '--r');
% %     box on; grid on; xlabel('\Delta knee_{joint accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); ylabel('rew [adim]'); hold off
% % 
% % subplot(2,2,4); plot(hipjerk, rewjointjerk,'b'); hold on; yline(0, '--g'); xline(0, '--r');
% %     box on; grid on; xlabel('\Delta hip_{joint accel} [rad/s^2]'); xlim([min(hipjerk), max(hipjerk)]); ylabel('rew [adim]'); hold off
% % 
% % fontsize(fig7, scale=1.2)   
% % 
% % %===================== 8. Reward Action Angular Jerk ======================%
% % fig8 = figure(8);   
% % sgtitle("Reward Action Angular Jerk")
% % subplot(2,2,1); surf(kneejerk, hipjerk, rewactjerk); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
% %     box on; grid on; xlabel('\Delta knee_{action accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action accel} [rad]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% % 
% % subplot(2,2,2); pcolor(kneejerk, hipjerk, rewactjerk);
% %     box on; grid on; xlabel('\Delta knee_{action accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action accel} [rad]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% % 
% % subplot(2,2,3); plot(kneejerk, rewactjerk, 'b'); hold on; yline(0, '--g'); xline(0, '--r');
% %     box on; grid on; xlabel('\Delta knee_{action accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); ylabel('rew [adim]'); hold off
% % 
% % subplot(2,2,4); plot(hipjerk, rewactjerk,'b'); hold on; yline(0, '--g'); xline(0, '--r');
% %     box on; grid on; xlabel('\Delta hip_{action accel} [rad/s^2]'); xlim([min(hipjerk), max(hipjerk)]); ylabel('rew [adim]'); hold off
% % 
% % fontsize(fig8, scale=1.2) 
% 
% %===================== 9. Reward Orientation =====================%
fig9 = figure(9);
sgtitle("Reward Orientation")
subplot(2,2,1); surf(prgrav_x, prgrav_x, reworient); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
    box on; grid on; xlabel('proj.grav_x [adim]'); xlim([min(prgrav_x) , max(prgrav_x)]); set(gca, 'XDir','reverse'); ylabel('proj.grav_y [adim]'); ylim([min(prgrav_y), max(prgrav_y)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off

subplot(2,2,2); pcolor(prgrav_x, prgrav_x, reworient);
   box on; grid on; xlabel('proj.grav_x [adim]'); xlim([min(prgrav_x) , max(prgrav_x)]); set(gca, 'XDir','reverse'); ylabel('proj.grav_y [adim]'); ylim([min(prgrav_y), max(prgrav_y)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;

subplot(2,2,3); plot(prgrav_x, reworient, 'b'); hold on; xline(0, '--r');
    box on; grid on; xlabel('proj.grav_x [adim]'); xlim([min(prgrav_x) , max(prgrav_x)]); ylabel('rew [adim]'); hold off

subplot(2,2,4); plot(prgrav_y, reworient,'b'); hold on; xline(0, '--r');
    box on; grid on; xlabel('proj.grav_y [adim]'); xlim([min(prgrav_y), max(prgrav_y)]); ylabel('rew [adim]'); hold off

fontsize(fig9, scale=1.2) 
% % 
% % %======================= 10. Reward Action Jerk Max ======================%
% % fig10 = figure(10); 
% % sgtitle("Reward Action Jerk Max")
% % subplot(2,2,1); surf(kneejerk, hipjerk, rewjerk_max); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
% %     box on; grid on; xlabel('\Delta knee_{action accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action accel} [rad/s^2]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% % 
% % subplot(2,2,2); pcolor(kneejerk, hipjerk, rewjerk_max);
% %     box on; grid on; xlabel('\Delta knee_{action accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action accel} [rad/s^2]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% % 
% % subplot(2,2,3); plot(kneejerk, rewjerk_max, 'b'); hold on; xline(0, '--r');
% %     box on; grid on; xlabel('\Delta knee_{action accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); ylabel('rew [adim]'); hold off
% % 
% % subplot(2,2,4); plot(hipjerk, rewjerk_max,'b'); hold on; xline(0, '--r');
% %     box on; grid on; xlabel('\Delta hip_{action accel} [rad/s^2]'); xlim([min(hipjerk), max(hipjerk)]); ylabel('rew [adim]'); hold off
% % 
% % fontsize(fig10, scale=1.2)   
% % 
% % %======================= 11. Reward Joint Jerk Max =======================%
% % fig11 = figure(11); 
% % sgtitle("Reward Joint Jerk Max")
% % subplot(2,2,1); surf(kneejerk, hipjerk, rewjjerk_max); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
% %     box on; grid on; xlabel('\Delta knee_{joint accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint accel} [rad/s^2]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% % 
% % subplot(2,2,2); pcolor(kneejerk, hipjerk, rewjjerk_max);
% %     box on; grid on; xlabel('\Delta knee_{joint accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint accel} [rad/s^2]'); ylim([min(hipjerk), max(hipjerk)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% % 
% % subplot(2,2,3); plot(kneejerk, rewjjerk_max, 'b'); hold on; xline(0, '--r');
% %     box on; grid on; xlabel('\Delta knee_{joint accel} [rad/s^2]'); xlim([min(kneejerk) , max(kneejerk)]); ylabel('rew [adim]'); hold off
% % 
% % subplot(2,2,4); plot(hipjerk, rewjjerk_max,'b'); hold on; xline(0, '--r');
% %     box on; grid on; xlabel('\Delta hip_{joint accel} [rad/s^2]'); xlim([min(hipjerk), max(hipjerk)]); ylabel('rew [adim]'); hold off
% % 
% % fontsize(fig11, scale=1.2)   
% 
% %======================= 12. Reward Action Accel Max ======================%
% fig12 = figure(12); 
% sgtitle("Reward Action Accel Max")
% subplot(2,2,1); surf(kneeaccel, hipaccel, rewaccel_max); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
%     box on; grid on; xlabel('\Delta knee_{action rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action rate} [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% 
% subplot(2,2,2); pcolor(kneeaccel, hipaccel, rewaccel_max);
%     box on; grid on; xlabel('\Delta knee_{action rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{action rate} [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% 
% subplot(2,2,3); plot(kneeaccel, rewaccel_max, 'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta knee_{action rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); ylabel('rew [adim]'); hold off
% 
% subplot(2,2,4); plot(hipaccel, rewaccel_max,'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta hip_{action rate} [rad/s]'); xlim([min(hipaccel), max(hipaccel)]); ylabel('rew [adim]'); hold off
% 
% fontsize(fig12, scale=1.2)   
% 
% %======================= 13. Reward Joint Accel Max =======================%
% fig13 = figure(13); 
% sgtitle("Reward Joint Accel Max")
% subplot(2,2,1); surf(kneeaccel, hipaccel, rewjaccel_max); hold on; planexy = surf(x, y, z,'FaceAlpha',0.5); set(planexy, 'Facecolor', "#4DBEEE", "EdgeColor", "none"); % Traccia il piano z = 0
%     box on; grid on; xlabel('\Delta knee_{joint rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint rate} [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); hold off
% 
% subplot(2,2,2); pcolor(kneeaccel, hipaccel, rewjaccel_max);
%     box on; grid on; xlabel('\Delta knee_{joint rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); set(gca, 'XDir','reverse'); ylabel('\Delta hip_{joint rate} [rad/s]'); ylim([min(hipaccel), max(hipaccel)]); set(gca, 'YDir','reverse'); zlabel('rew [adim]'); axis equal;
% 
% subplot(2,2,3); plot(kneeaccel, rewjaccel_max, 'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta knee_{joint rate} [rad/s]'); xlim([min(kneeaccel) , max(kneeaccel)]); ylabel('rew [adim]'); hold off
% 
% subplot(2,2,4); plot(hipaccel, rewjaccel_max,'b'); hold on; xline(0, '--r');
%     box on; grid on; xlabel('\Delta hip_{joint rate} [rad/s]'); xlim([min(hipaccel), max(hipaccel)]); ylabel('rew [adim]'); hold off
% 
% fontsize(fig13, scale=1.2)   

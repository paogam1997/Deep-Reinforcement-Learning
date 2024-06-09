% File analogo di Domrandom_test in cui tutti i conti dei results e print
% sono ottimizzati in cicli e non fatti manualmente
clc;
% Va aggiunta la path la cartella con le funzioni di supporto create
addpath ('Funzioni Utili')

%% PARAMETRI GLOBALI
% Definizione di grandezze di simulazione e task
policy_step_dt = 0.01;               % [s] durata di uno step della policy
simulation_dt = 0.005;               % dentro un ciclo della policy si effettua una simulazione di 5ms dopo la pre_physics e un'altro dopo la post_physics
episode_lenght = 4000;               % [policy_step] = 40 s (effettivi)(e non 20 s come scritto in yaml) (questa differenza non fa cambiare i ragionamenti successivi)
 base_z_initial_condition = 0.35;    % [m] altezza iniziale a cui è fatto spawnare Mulinex
%base_z_initial_condition = 0.365;   % con questa non va subito in collisione
toll_check_collision_orient = 0.01;  % [m] tolleranza minima da lasciare per non rilevare contatto con terreno nel check per la orientazione
massa_Mulinex = 5.5;                 % [kg] massa del robot
g = 9.81;                            % [m/s^2] accelerazione di gravità

%% INIZIALIZZAZIONE
% Inizializzazione della struttura dei parametri da randomizzare
Params2Rand = [];
Params2Rand.actions = [];              Params2Rand.gravity = [];
Params2Rand.stat_terr_friction = [];   Params2Rand.dyn_terr_friction = [];
Params2Rand.restitution = [];
Params2Rand.orientation = [];          Params2Rand.stiffness = [];
Params2Rand.damping = [];              Params2Rand.joint_friction = [];
Params2Rand.lower_dof_limits = [];     Params2Rand.upper_dof_limits = [];
Params2Rand.joint_armatures = [];      Params2Rand.max_efforts = [];
Params2Rand.joint_efforts = [];
% <== Aggiungere qua eventuali nuovi parametri

% Inizializzazione strutture associate ai sigoli parametri
actions = [];               gravity = [];
stat_terr_friction = [];    dyn_terr_friction = [];
restitution = [];
orientation = [];           stiffness = [];
damping = [];               joint_friction = [];
lower_dof_limits = [];      upper_dof_limits = [];
joint_armatures = [];       max_efforts = [];
joint_efforts = [];
% <== Aggiungere qua eventuali nuovi parametri

% Creazione di una lista dei parametri da randomizzare
params = fieldnames(Params2Rand);
% Creazione delle prime chiavi del dizionario associato ai risultati
index_results = ["Freq [s]", "Freq [p_steps]"];

%% DEFINIZIONI RANDOMIZZAZIONI
% Ogni riga di dr_info corrisponde ad una riga di distribution_parameters
% (ci metto solo le varianze delle gaussiane)
% ed equivale ad 1 metodo di randomizzazione
%-------------------------------- ACTIONS --------------------------------%
actions.dr_info = [["on_reset", "additive", "gaussian"];
                   ["on_interval", "additive", "gaussian"]];

actions.distribution_parameters = [DeltaToVar_angoli(0); 
                                   DeltaToVar_angoli(0.314)];

actions.frequency_interval = SecondToPolicyStepDuration(2.5, policy_step_dt);
actions.type = "angle";     actions.unity = "[deg]";
%-------------------------------- GRAVITY ---------------------------------%
gravity.dr_info = [["on_reset", "additive", "uniform"];
                   ["on_interval", "None", "None"]];

gravity.distribution_parameters = [[0, 0.375*g]; % può solo incrementare
                                   [0,0]]; % Secondo valore è per eventuale on_interval
gravity.frequency_interval = 0;
gravity.type = "general";   gravity.unity = "[m/s^2]";
    % INFO: -Mulinex pesa 5.5-6 kg
    %       - regge al massimo altri 5-10 kg
    %       - braccio accessorio pesa circa 1-1.5 kg
%------------------------ STATIC TERRAIN FRICTION ------------------------%
stat_terr_friction.dr_info = [["on_reset", "direct", "uniform"];
                              ["on_interval", "scaling", "uniform"]];
stat_terr_friction.distribution_parameters = [[0.7, 1];
                                              [0.45, 1.3]];
stat_terr_friction.frequency_interval = 2350;
stat_terr_friction.type = "general";    stat_terr_friction.unity = "[adim]";
    % INFO: - è la prima componente del vettore unico "material_property"
    %         usato nello yaml.
%----------------------- DYNAMIC TERRAIN FRICTION ------------------------%
dyn_terr_friction.dr_info = [["on_reset", "direct", "uniform"];
                              ["on_interval", "scaling", "uniform"]];
dyn_terr_friction.distribution_parameters = [[0.7, 1];
                                             [0.21, 1.3]];
dyn_terr_friction.frequency_interval = 2350;
dyn_terr_friction.type = "general";    dyn_terr_friction.unity = "[adim]";
    % INFO: - è la seconda componente del vettore unico "material_property"
    %         usato nello yaml.
%--------------------- TERRAIN ELASTIC RESTITUTION -----------------------%
restitution.dr_info = [["on_reset", "direct", "uniform"];
                       ["on_interval", "scaling", "uniform"]];
restitution.distribution_parameters = [[0.5, 1];
                                       [0.5, 1.5]];
restitution.frequency_interval = 2350;
restitution.type = "general";    restitution.unity = "[adim]";
    % INFO: - è la terza componente del vettore unico "material_property"
    %         usato nello yaml.
    %       - coeff. di restituzione elastica è il raporto tra la velocità
    %         relativa dopo l'impatto e la velocità relativa prima dell'impatto.
    %       - 0<e<1, e=1 urto totalmente elastico --> terreno perfettamente
    %         rigido, e=0 urto totalmente anaelastico
%------------------------------ ORIENTATION ------------------------------%
% Versione base
% orientation.dr_info = [["on_reset", "additive", "gaussian"];
%                        ["on_interval", "additive", "gaussian"]];

%orientation.distribution_parameters = [[0.00338, 0.00338, 0.01014]; 
%                                       [0.0000001322, 0.0000001542, 0.0000000562]];
    % Versione prova
    orientation.dr_info = [["on_reset", "additive", "uniform"]; 
                           ["on_interval", "additive", "uniform"]];

    % Per additive uniform il parametro da dare è il max dell'intervallo
    % [min, max], supposto simmetrico
    orientation.distribution_parameters = [[deg2radiant(0.5), deg2radiant(0.5), deg2radiant(0.10)]; 
                                           [deg2radiant(0.0), deg2radiant(0.0), deg2radiant(0.0)]]; % sono come roll pitch e yaw ma yaml vuole angoli di Eulero
orientation.frequency_interval = 100;
orientation.type = "angle";     orientation.unity = "[deg]";
orientation.form = "RPY";
    % INFO: - Da .yaml non funziona, implementata solo on_reset direttamente
    %         da codice in .py
%------------------------------- STIFFNESS -------------------------------%
stiffness.dr_info = [["on_reset", "scaling", "uniform"];
                     ["on_interval", "scaling", "uniform"]];

stiffness.distribution_parameters = [[0.9, 1.1];
                                     [0.8, 1.2]];
stiffness.frequency_interval = 2000;
stiffness.type = "general";     stiffness.unity = "[Nm/rad]";
    % INFO: - Valore nominale 45.5 Nm/rad
%------------------------------- DAMPING ---------------------------------%
damping.dr_info = [["on_reset", "scaling", "uniform"];
                   ["on_interval", "scaling", "uniform"]];

damping.distribution_parameters = [[0.8, 1.2]; 
                                   [0.9, 1.1]];
damping.frequency_interval = 1750;
damping.type = "general";       damping.unity = "[Nms^2/rad]";
    % INFO: - Valore nominale 0.26 Nms^2/rad
%----------------------------- JOINT FRICTION ----------------------------%
joint_friction.dr_info = [["on_reset", "scaling", "uniform"]; 
                          ["on_interval", "scaling", "uniform"]];

joint_friction.distribution_parameters = [[0.85, 1.15]; 
                                          [0.82, 1.18]];
joint_friction.frequency_interval = 1300;
joint_friction.type = "general";    joint_friction.unity = "[adim]";
%------------------------------ LOWER DOF LIMITS -------------------------%
lower_dof_limits.dr_info = [["on_reset", "additive", "gaussian"]; 
                            ["on_interval", "additive", "gaussian"]];
% Versione base
%lower_dof_limits.distribution_parameters = [0.00338; 
%                                           0.000846];
    % Versione prova
    lower_dof_limits.distribution_parameters = [DeltaToVar_angoli(10);
                                                DeltaToVar_angoli(6)];
lower_dof_limits.frequency_interval = 1850;
lower_dof_limits.type = "angle";    lower_dof_limits.unity = "[deg]";
%---------------------------- UPPER DOF LIMITS ----------------------------%
upper_dof_limits.dr_info = [["on_reset", "additive", "gaussian"]; 
                            ["on_interval", "additive", "gaussian"]];
% Versione base
% upper_dof_limits.distribution_parameters = [0.00338; 
%                                             0.000846];
    % Versione prova
    upper_dof_limits.distribution_parameters = [DeltaToVar_angoli(10); 
                                                DeltaToVar_angoli(6)];
upper_dof_limits.frequency_interval = 1850;
upper_dof_limits.type = "angle";    upper_dof_limits.unity = "[deg]";
%--------------------------- JOINT ARMATURES -----------------------------%
joint_armatures.dr_info = [["on_reset", "direct", "unifom"]; 
                           ["on_interval", "None", "None"]];

joint_armatures.distribution_parameters = [[4*10^-3, 4*10^-3]; 
                                           [0, 0]];
joint_armatures.frequency_interval = 950;
joint_armatures.type = "general";   joint_armatures.unity = "[kg*m^2]";
    % INFO: - Valore nominale 2*10^-3 kgm^2
%------------------------------ MAX EFFORTS ------------------------------%
max_efforts.dr_info = [["on_reset", "additive", "uniform"]; 
                       ["on_interval", "None", "None"]];
                         
max_efforts.distribution_parameters = [[-1, 0];
                                         [0, 0]];
max_efforts.frequency_interval = 0;
max_efforts.type = "general"; max_efforts.unity = "[Nm]";
%----------------------------- JOINT EFFORTS -----------------------------%
joint_efforts.dr_info = [["on_reset", "scaling", "uniform"]; 
                         ["on_interval", "scaling", "uniform"]];

joint_efforts.distribution_parameters = [[1,1];
                                         [0.85, 1.15]];
joint_efforts.frequency_interval = 2600;
joint_efforts.type = "general"; joint_efforts.unity = "[Nm]";

% <== Aggiungere qua eventuali nuovi parametri

% Associo le strutture dei parametri ai campi della macrostruttura
% Params2Rand
Params2Rand.actions = actions;
Params2Rand.gravity = gravity;
Params2Rand.stat_terr_friction = stat_terr_friction;
Params2Rand.dyn_terr_friction = dyn_terr_friction;
Params2Rand.restitution = restitution;
Params2Rand.orientation = orientation;
Params2Rand.stiffness = stiffness;
Params2Rand.damping = damping;
Params2Rand.joint_friction = joint_friction;
Params2Rand.lower_dof_limits = lower_dof_limits;
Params2Rand.upper_dof_limits = upper_dof_limits;
Params2Rand.joint_armatures = joint_armatures;
Params2Rand.max_efforts = max_efforts;
Params2Rand.joint_efforts = joint_efforts;
% <== Aggiungere qua eventuali nuovi parametri

%% ======== CALCOLO DEI VARI TERMINI UTILI, FINIRANNO IN .RESULTS ====== %%
% Da qua in poi non toccare più nulla, i calcoli sono tutti automatizzati e
% si adattano da soli alle grandezze che inserisci. Si modifica qua sotto
% solo se si introducono grandezze vettoriali come orientation, in tal caso
% seguire quanto fatto con lei.

%% CONVERSIONE FREQUENCY_INTERVAL IN S

freq_s = zeros(1, numel(params));   % Inizializzazione variabile
for i=1:numel(params)               % Calcolo la variabile per ogni elemento di params = Params2Rand.(params{i})
    freq_s(i)= PolicyStepDurationToSeconds(Params2Rand.(params{i}).frequency_interval, policy_step_dt);                    % Conversione della frequenza data da policy steps in secondi
    Params2Rand.(params{i}).results = dictionary(index_results, [freq_s(i), Params2Rand.(params{i}).frequency_interval]);  % Creo il dizionario e lo inizializzo con le frequenze in secondi e policy steps
end 

%% CALCOLO MAX VARIATION ON_RESET

for i=1:numel(params)  % Calcolo la variabile per ogni elemento di params = Params2Rand.(params{i})
    if ~ strcmp(params{i}, 'orientation')  % Eseguo il calcolo a meno che il parametro non sia "orientation" dato che è vettoriale
        Params2Rand.(params{i}).results("MaxD_reset") = MaxDelta(Params2Rand.(params{i}), "on_reset");   % Calcolo della massima variazione causata dalla randomizzazione on_reset
        if strcmp(params{i}, 'gravity')   % Se il parametro è la gravità creo due voci nel dizionario in più per riportare la stessa variazione in [g] e in [kg] riferendosi al peso di Mulinex
            Params2Rand.(params{i}).results("MaxD_reset") = Params2Rand.(params{i}).results("MaxD_reset")/9.81;                % Definizione della variazione come multiplo di g
            Params2Rand.(params{i}).results("MaxD_reset [kg]") = Params2Rand.(params{i}).results("MaxD_reset")*massa_Mulinex;  % Definizione della variazione come multiplo della massa di Mulinex
        end
    end
end
    
%------- I results per orientation vanno comunque calcolati a mano -------%
orientation_MaxVariationReset = zeros(1, 3);                                         % Inizializzazione della variabile
if orientation.dr_info(1,2) == "additive" && orientation.dr_info(1,3) == "gaussian"  % Nel caso si scelga una randomizzazione additiva con ddp gaussiana
    for j = 1:size(orientation.distribution_parameters, 2)                           % Iterazione su roll, pitch e yaw
        orientation_MaxVariationReset (1, j) = VarToDelta_angoli(orientation.distribution_parameters(1,j));  % Calcolo della massima variazione on_reset
    end
elseif orientation.dr_info(1,2) == "additive" && orientation.dr_info(1,3) == "uniform"
    for j = 1:size(orientation.distribution_parameters, 2)
        orientation_MaxVariationReset (1, j) = orientation.distribution_parameters(1,j) * (180/pi);
    end
end
Params2Rand.orientation.results("roll_MaxD_reset") = orientation_MaxVariationReset(1,1);
Params2Rand.orientation.results("pitch_MaxD_reset") = orientation_MaxVariationReset(1,2);
Params2Rand.orientation.results("yaw_MaxD_reset") = orientation_MaxVariationReset(1,3);

%% FINAL DELTA E REPETITION

% Action e gravity non hanno questi parametri --> for salta i primi 2
% elementi di params
for i=1:numel(params)
    if ~ strcmp(params{i}, 'orientation')
        [Params2Rand.(params{i}).results("FinalDelta"), ~, Params2Rand.(params{i}).results("Repetition")]  = FinalMaxDelta(Params2Rand.(params{i}), episode_lenght);
        if strcmp(params{i}, 'gravity')
            Params2Rand.(params{i}).results("FinalDelta") = Params2Rand.(params{i}).results("FinalDelta")/9.81; % Misuro la variazione in g
            Params2Rand.(params{i}).results("FinalDelta [kg]") = Params2Rand.(params{i}).results("FinalDelta")*massa_Mulinex;
        end
        if Params2Rand.(params{i}).dr_info(1, 2) == "additive" && Params2Rand.(params{i}).dr_info(1, 3) == "gaussian"
            Params2Rand.(params{i}).results("OnInt Incr") = Params2Rand.(params{i}).results("FinalDelta") - Params2Rand.(params{i}).results("MaxD_reset");
        elseif Params2Rand.(params{i}).dr_info(1, 2) == "scaling" && Params2Rand.(params{i}).dr_info(1, 3) == "uniform"
            Params2Rand.(params{i}).results("OnInt Incr") = Params2Rand.(params{i}).results("FinalDelta") / Params2Rand.(params{i}).results("MaxD_reset");
        end
        
            
    end
end

%------ Orientation riciede calcolo separato perchè ha 3 componenti ------%
% I results per orientation vanno comunque calcolati a mano
% Params2Rand.orientation.results("FinalDelta") = 0.0; % non so come cancellare la voce quindi la setto a 0
% Params2Rand.orientation.results("OnInt Incr") = 0.0;

int_orient = floor(episode_lenght /orientation.frequency_interval);
delta_orient = zeros(int_orient+1, 3);
for j = 1:size(orientation.distribution_parameters, 2) % Itero su roll, pitch e yaw
    if orientation.dr_info(1,2) == "additive" && orientation.dr_info(1,3) == "gaussian"
        % Ogni angolo ha una colonna
        % Prima riga è il delta dato da on_reset
        delta_orient(1,j) = VarToDelta_angoli(orientation.distribution_parameters(1,j));
        % Riempio altre righe con i delta di on_interval
        for i = 1:int_orient
            delta_orient(i+1,j) = VarToDelta_angoli(orientation.distribution_parameters(2,j));
        end
    elseif orientation.dr_info(1,2) == "additive" && orientation.dr_info(1,3) == "uniform"
        delta_orient(1,j) =  orientation.distribution_parameters(1,j) * (180/pi);
        delta_orient(2:end,j) = orientation.distribution_parameters(2,j) * (180/pi);
    end
end
if orientation.dr_info(1, 2) == "additive" && orientation.dr_info(1, 3) == "gaussian"
    complete_delta_orient = sum(delta_orient, 1); % vettore 1 riga (risultati) e 3 colonne (roll, pitch, yaw)
elseif orientation.dr_info(1, 2) == "scaling" && orientation.dr_info(1, 3) == "uniform"
    complete_delta_orient = prod(delta_orient, 1); % vettore 1 riga (risultati) e 3 colonne (roll, pitch, yaw)
elseif orientation.dr_info(1,2) == "additive" && orientation.dr_info(1,3) == "uniform"
    complete_delta_orient = sum(delta_orient, 1); % vettore 1 riga (risultati) e 3 colonne (roll, pitch, yaw)
end
delta_orient_ = delta_orient(:,:);
repetition_orient = int_orient;

Params2Rand.orientation.results("Repetition") = repetition_orient;
Params2Rand.orientation.results("roll_FinalDelta") = complete_delta_orient(1,1);
Params2Rand.orientation.results("pitch_FinalDelta") = complete_delta_orient(1,2);
Params2Rand.orientation.results("yaw_FinalDelta") = complete_delta_orient(1,3);
if Params2Rand.orientation.dr_info(1, 2) == "additive" && Params2Rand.orientation.dr_info(1, 3) == "gaussian"
    Params2Rand.orientation.results("roll OnInt Incr") = Params2Rand.orientation.results("roll_FinalDelta") - Params2Rand.orientation.results("roll_MaxD_reset");
    Params2Rand.orientation.results("pitch OnInt Incr") = Params2Rand.orientation.results("pitch_FinalDelta") - Params2Rand.orientation.results("pitch_MaxD_reset");
    Params2Rand.orientation.results("yaw OnInt Incr") = Params2Rand.orientation.results("yaw_FinalDelta") - Params2Rand.orientation.results("yaw_MaxD_reset");
elseif Params2Rand.orientation.dr_info(1, 2) == "scaling" && Params2Rand.orientation.dr_info(1, 3) == "uniform"
    Params2Rand.orientation.results("roll OnInt Incr") = Params2Rand.orientation.results("roll_FinalDelta") / Params2Rand.orientation.results("roll_MaxD_reset");
    Params2Rand.orientation.results("pitch OnInt Incr") = Params2Rand.orientation.results("pitch_FinalDelta") / Params2Rand.orientation.results("pitch_MaxD_reset");
    Params2Rand.orientation.results("yaw OnInt Incr") = Params2Rand.orientation.results("yaw_FinalDelta") / Params2Rand.orientation.results("yaw_MaxD_reset");
elseif Params2Rand.orientation.dr_info(1, 2) == "additive" && Params2Rand.orientation.dr_info(1, 3) == "uniform"
    Params2Rand.orientation.results("roll OnInt Incr") = Params2Rand.orientation.results("roll_FinalDelta") - Params2Rand.orientation.results("roll_MaxD_reset");
    Params2Rand.orientation.results("pitch OnInt Incr") = Params2Rand.orientation.results("pitch_FinalDelta") - Params2Rand.orientation.results("pitch_MaxD_reset");
    Params2Rand.orientation.results("yaw OnInt Incr") = Params2Rand.orientation.results("yaw_FinalDelta") - Params2Rand.orientation.results("yaw_MaxD_reset");
end

%% PRINT DEI RESULTS

for i=1:numel(params)
    disp("======================================================")
    fprintf("---------------------- %s -----------------------\n", params{i});
    disp("======================================================")
    fprintf("Measurament unity      %s\n", Params2Rand.(params{i}).unity)
    if Params2Rand.(params{i}).dr_info(1,2) == "None"
        fprintf("Type of randomization:     %s, %s\n", Params2Rand.(params{i}).dr_info(2,2), Params2Rand.(params{i}).dr_info(2,3))
        fprintf("Only ''%s'' randomization\n", Params2Rand.(params{i}).dr_info(2,1))
    elseif Params2Rand.(params{i}).dr_info(2,2) == "None"
        fprintf("Type of randomization:     %s, %s\n", Params2Rand.(params{i}).dr_info(1,2), Params2Rand.(params{i}).dr_info(1,3))
        fprintf("Only ''%s'' randomization\n", Params2Rand.(params{i}).dr_info(1,1))
    else
        fprintf("On_reset randomization:     %s, %s\n", Params2Rand.(params{i}).dr_info(1,2), Params2Rand.(params{i}).dr_info(1,3))
        fprintf("On_interval randomization:     %s, %s\n", Params2Rand.(params{i}).dr_info(2,2), Params2Rand.(params{i}).dr_info(2,3))
    end
    Params2Rand.(params{i}).results
    if strcmp(params{i}, 'gravity')
        disp("Le variazioni sono in [g = gravità]")
    end
    if strcmp(params{i}, 'orientation')
        if orientation.dr_info(1,3) == "uniform"
            euler = rpy2eul(rad2deg(orientation.distribution_parameters(1,1)), rad2deg(orientation.distribution_parameters(1,2)), rad2deg(orientation.distribution_parameters(1,3)));
            euler_int = rpy2eul(rad2deg(orientation.distribution_parameters(2,1)), rad2deg(orientation.distribution_parameters(2,2)), rad2deg(orientation.distribution_parameters(2,3)));
        elseif orientation.dr_info(1,3) == "gaussian"
            euler = rpy2eul(VarToDelta_angoli(orientation.distribution_parameters(1,1)), VarToDelta_angoli(orientation.distribution_parameters(1,2)), VarToDelta_angoli(orientation.distribution_parameters(1,3)));
            euler_int = rpy2eul(VarToDelta_angoli(orientation.distribution_parameters(2,1)), VarToDelta_angoli(orientation.distribution_parameters(2,2)), VarToDelta_angoli(orientation.distribution_parameters(2,3)));
        end
        if orientation.form == "RPY"
            disp("Max Delta in Eulero:")
            if orientation.distribution_parameters(1,1) ~= 0 && orientation.distribution_parameters(1,2) ~= 0
                fprintf("   on_reset: [%s, %s, %s]\n", [euler(:)])
            else
                fprintf("   on_reset: n.d.\n")
            end
            if orientation.distribution_parameters(2,1) ~= 0 && orientation.distribution_parameters(2,2) ~= 0
                fprintf("   on_interval: [%s, %s, %s]\n", [euler_int(:)])
            else
                fprintf("   on_interval: n.d.\n")
            end
        end
        if orientation.dr_info(1,2) == "additive" && orientation.dr_info(1,3) == "uniform"
            if orientation.distribution_parameters(1,1) ~= 0 && orientation.distribution_parameters(1,2) ~= 0
                min_quat_reset = eul2quat(-orientation.distribution_parameters(1,:));
                max_quat_reset = eul2quat(orientation.distribution_parameters(1,:));
                fprintf("Quaternione minimo in on_reset: [%s, %s, %s, %s]\n", [min_quat_reset(:)])
                fprintf("Quaternione massimo in on_reset: [%s, %s, %s, %s]\n", [max_quat_reset(:)])
            end
        end
        if orientation.dr_info(1,2) == "additive" && orientation.dr_info(1,3) == "uniform"
            if orientation.distribution_parameters(2,1) ~= 0 && orientation.distribution_parameters(2,2) ~= 0
                min_quat_int = eul2quat(-orientation.distribution_parameters(2,:));
                max_quat_int = eul2quat(orientation.distribution_parameters(2,:));
                fprintf("Quaternione minimo in on_interval: [%s, %s, %s, %s]\n", [min_quat_int(:)])
                fprintf("Quaternione massimo in on_interval: [%s, %s, %s, %s]\n", [max_quat_int(:)])
            end
        end

        % Ingressi di CheckCollision.. sono in gradi
        info_final = CheckCollisionRandOrient(Params2Rand.orientation.results("roll_FinalDelta"), Params2Rand.orientation.results("pitch_FinalDelta"), Params2Rand.orientation.results("yaw_FinalDelta"), base_z_initial_condition, toll_check_collision_orient);
        info_on_reset = CheckCollisionRandOrient(Params2Rand.orientation.results("roll_MaxD_reset"), Params2Rand.orientation.results("pitch_MaxD_reset"), Params2Rand.orientation.results("yaw_MaxD_reset"), base_z_initial_condition, toll_check_collision_orient);
        disp(info_on_reset)
        disp(info_final)
        disp("  ")
    end
end

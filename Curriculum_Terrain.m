max_level = 10;                    % massimo livello di difficoltà
act_level = 0:max_level;           % livello attuale
difficulty = act_level/max_level;  % parametro di diffoltà che aumenta ad ogni livello, usato per rendere più ostici gli ostacoli

%% Tipi di terreni
    % Base
slope_b = difficulty * 0.4;
step_height_b = 0.05 + 0.175 * difficulty;
discrete_obstacles_height_b = 0.025 + difficulty * 0.15;
stepping_stones_size_b = 2 - 1.8 * difficulty;
wave_amplitude = 0.04 + 0.12 * difficulty;
    
% Double level
max_level1 = 20;   act_level1 = 0:max_level1;  difficulty1 = act_level1/max_level1;  
slope_1 = difficulty1 * 0.4;
step_height_1 = 0.05 + 0.175 * difficulty1;
discrete_obstacles_height_1 = 0.025 + difficulty1 * 0.15;
stepping_stones_size_1 = 2 - 1.8 * difficulty1;
wave_amplitude_1 = 0.04 + 0.12 * difficulty1;
    
% Modified terrain (Easy Terrain versione vecchia)
max_level2 = 30;   act_level2 = 0:max_level2;  difficulty2 = act_level2/max_level2;  
slope_2 = difficulty2 * 0.4;
step_height_2 = 0.025 + 0.175 * difficulty2;
discrete_obstacles_height_2 = 0.0125 + difficulty2 * 0.15;
stepping_stones_size_2 = 2 - 1.8 * difficulty2;
wave_amplitude_2 = 0.04 + 0.12 * difficulty;
    
% Extra Easy Terrain 3
max_level3 = 20;   act_level3 = 0:max_level3;  difficulty3 = act_level3/max_level3;  
slope_3 = difficulty3 * 0.388;
step_height_3 = 0.01 + 1/18 * difficulty3;
discrete_obstacles_height_3 = 0.005 + difficulty3 * 1/30;
stepping_stones_size_3 = 2 - 1.8 * difficulty3;
wave_amplitude_3 = 0.01 + 1/39 * difficulty3;

% Easy Terrain 2
max_level4 = 20;   act_level4 = 0:max_level4;  difficulty4 = act_level4/max_level4;  
slope_4 = difficulty4 * 0.4;
step_height_4 = 0.025 + 0.175 * difficulty4;
discrete_obstacles_height_4 = 0.0125 + difficulty4 * 0.15;
stepping_stones_size_4 = 2 - 1.8 * difficulty4;
wave_amplitude_4 = 0.04 + 0.12 * difficulty4;

% Easy Terrain 3
max_level5 = 20;   act_level5 = 0:max_level5;  difficulty5 = act_level5/max_level5;  
slope_5 = difficulty5 * 0.4;
step_height_5 = 0.01 + 0.2 * difficulty5;
discrete_obstacles_height_5 = 0.005 + difficulty5 * 0.1;
stepping_stones_size_5 = 2 - 1.8 * difficulty5;
wave_amplitude_5 = 0.01 + 0.1 * difficulty5;

%% Plot
fig1 = figure(1); sgtitle("Terrain Parameters")
subplot(2,2,1); 
 %   plot(act_level, slope_b*180/pi, 'bo-'); hold on; 
 %   plot(act_level1, slope_1*180/pi, 'ro-');
 %   plot(act_level2, slope_2*180/pi, 'go-');
    plot(act_level3, slope_3*180/pi, 'x-', LineWidth=2); hold on;
    plot(act_level4, slope_4*180/pi, '*-', LineWidth=2);
    plot(act_level5, slope_5*180/pi, 'o-', LineWidth=2);
    box on; grid on; title("Slopes"); xlabel("Actual Level [adim]"); ylabel("[deg]"); legend(["EET3", "ET2","ET3"],"Location","southeast"); hold off;
    set(gca, 'FontSize', 16);
subplot(2,2,2); 
 %   plot(act_level, step_height_b*100, 'bo-'); hold on; 
 %   plot(act_level1, step_height_1*100, 'ro-');
  %  plot(act_level2, step_height_2*100, 'go-');
    plot(act_level3, step_height_3*100, 'x-', LineWidth=2); hold on;
    plot(act_level4, step_height_4*100, '*-', LineWidth=2);
    plot(act_level5, step_height_5*100, 'o-', LineWidth=2);
    box on; grid on; title("Steps height"); xlabel("Actual Level [adim]"); ylabel("[cm]"); legend(["EET3", "ET2","ET3"],"Location","northwest"); hold off;
    set(gca, 'FontSize', 16);
subplot(2,2,3); 
%    plot(act_level, discrete_obstacles_height_b*100, 'bo-'); hold on;
%    plot(act_level1, discrete_obstacles_height_1*100, 'ro-');
%    plot(act_level2, discrete_obstacles_height_2*100, 'go-');
    plot(act_level3, discrete_obstacles_height_3*100, 'x-', LineWidth=2); hold on;
    plot(act_level4, discrete_obstacles_height_4*100, '*-', LineWidth=2);
    plot(act_level5, discrete_obstacles_height_5*100, 'o-', LineWidth=2);
    box on; grid on; title("Discrete obstacles height"); xlabel("Actual Level [adim]"); ylabel("[cm]"); legend(["EET3", "ET2","ET3"],"Location","northwest"); hold off;
    set(gca, 'FontSize', 16);
subplot(2,2,4); 
 %   plot(act_level, stepping_stones_size_b, 'bo-'); hold on;
%    plot(act_level1, stepping_stones_size_1, 'ro-');
  %  plot(act_level2, stepping_stones_size_2, 'go-');
    plot(act_level3, stepping_stones_size_3, 'x-', LineWidth=2); hold on;
    plot(act_level4, stepping_stones_size_4, '*-', LineWidth=2);
    plot(act_level5, stepping_stones_size_5, 'o-', LineWidth=2);
    box on; grid on; title("Larghezza pedane"); xlabel("livello attuale [adim]"); ylabel("stepping stones size [m]"); legend(["EET3",  "ET2","ET3"]); hold off;
    set(gca, 'FontSize', 16);
subplot(2, 2, 4);
    plot(act_level3, wave_amplitude_3*100, 'x-', LineWidth=2); hold on;
    plot(act_level4, wave_amplitude_4*100, '*-', LineWidth=2);
    plot(act_level5, wave_amplitude_5*100, 'o-', LineWidth=2);
    box on; grid on; title("Waves Aplitude"); xlabel("Actual Level [adim]"); ylabel("[cm]"); legend(["EET3", "ET2", "ET3"],"Location","northwest"); hold off;
    set(gca, 'FontSize', 16);
set(gca, 'FontSize', 16);
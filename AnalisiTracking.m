% ==================== Inserimento dati sperimentali ==================== %
Policy_name = {'RGR2', 'ET2x2Res', 'ET3x2Res', 'DR8EET3x2Res', 'BGR4ET3', 'EGREET3'}; % Policy testate
color = ["#0072BD", "#D95319", "#77AC30", "#7E2F8E", "#A2142F", "#EDB120"];                      % Associo ad ogni policy un colore per distinguerla nei grafici
vx_com = [0.3, 0.5, 0.7, 1, -0.3, -0.5, -0.7, -1];                         % [m/s] velocità lineari su x da testare
index_test_results = {'test_1', 'test_2', 'test_3'};   % Nomi degli esperimenti fatti per ogni velocità lineare
spost_des = 2;                                                             % [m] spostamento lineare lungo x di riferimento
% Ogni riga corrisponde ad un valore di velocità, ogni colonna ad un test
Results.RGR2.spost = [[0.65, 0.65, 0.74];           % vx=0.3
                       [1.168, 1.168, 1.205];         % vx=0.5
                       [1.675, 1.68, 1.72];     % vx=0.7
                       [1.54, 1.58, 1.41];     % vx=1
                       [0.65, 0.67, 0.78];     % vx=-0.3
                       [1.07, 0.855, 0.815];     % vx=-0.5
                       [0.87, 1.04, 1.155];     % vx=-0.7
                       [1.09, 1.17, 1.335]      % vx=-1
                       ];                    % [m] Spostamento longitudinale ottenuto
Results.RGR2.dev_or = [[0, 0.105, 0.105];       % vx=0.3
                       [0.25, 0.28, 0.185];      % vx=0.5
                       [0.04, 0.1, 0.09];  % vx=0.7
                       [0.04, 0.085, 0.07];  % vx=1
                       [0, 0.09, 0.095];  % vx=-0.3
                       [0.03, 0.03, 0.105];  % vx=-0.5
                       [0.063, 0.1, 0.035];  % vx=-0.7
                       [0.02, 0.05, 0]   % vx=-1
                       ];                           % [m] Deviazione laterale ottenuta

Results.ET2x2Res.spost = [[1.31, 1.385, 1.45];          % vx=0.3
                          [1.523, 1.48, 1.74];         % vx=0.5
                          [1.40, 1.40, 1.35];     % vx=0.7
                          [1.288, 1.323, 1.34];     % vx=1
                          [0.35, 0.444, 0.395];     % vx=-0.3
                          [1.14, 1.12, 1.11];     % vx=-0.5
                          [1.425, 1.50, 1.53];     % vx=-0.7
                          [0, 0, 0]      % vx=-1
                          ];                 % [m] Spostamento longitudinale ottenuto
Results.ET2x2Res.dev_or = [[0.43, 0.19, 0.195];       % vx=0.3
                           [0.05, 0.094, 0];      % vx=0.5
                           [0.03, 0.07, 0];  % vx=0.7
                           [0.095, 0.115, 0.09];  % vx=1
                           [0.058, 0.09, 0.15];  % vx=-0.3
                           [0.812, 1.01, 1.24];  % vx=-0.5
                           [0.43, 0.58, 0.675];  % vx=-0.7
                           [0, 0, 0]   % vx=-1
                           ];                        % [m] Deviazione laterale ottenuta

Results.ET3x2Res.spost = [[0.93, 1.30, 2.10];          % vx=0.3
                          [2.57, 1.74, 1.52];         % vx=0.5
                          [1.385, 1.635, 1.025];     % vx=0.7
                          [1.81, 1.74, 1.62];     % vx=1
                          [1.18, 1.205, 1.24];     % vx=-0.3
                          [0.82, 0.9, 1];     % vx=-0.5
                          [1.265, 1.345, 1.45];     % vx=-0.7
                          [1.38, 1.3, 1.34]      % vx=-1
                          ];                % [m] Spostamento longitudinale ottenuto

Results.ET3x2Res.dev_or = [[0.295, 0.6, 0.755];       % vx=0.3
                           [0.3, 0.26, 0.17];      % vx=0.5
                           [0.29, 0.035, 0.285];  % vx=0.7
                           [0, 0.043, 0.24];  % vx=1
                           [0.22, 0.27, 0.23];  % vx=-0.3
                           [0.08, 0.02, 0.08];  % vx=-0.5
                           [0.315, 0.305, 0.41];  % vx=-0.7
                           [0.36, 0.23, 0.13]   % vx=-1
                           ];                        % [m] Deviazione laterale ottenuta

Results.DR8EET3x2Res.spost = [[1, 1.055, 1.235];          % vx=0.3
                          [1.275, 1.49, 2.14];         % vx=0.5
                          [1.83, 1.75, 1.85];     % vx=0.7
                          [1.76, 1.79, 1.8];     % vx=1
                          [1.22, 1.24, 1.27];     % vx=-0.3
                          [1.19, 1.20, 1.26];     % vx=-0.5
                          [1.05, 1.105, 1.135];     % vx=-0.7
                          [1.06, 1.09, 1.21]      % vx=-1
                          ];                % [m] Spostamento longitudinale ottenuto

Results.DR8EET3x2Res.dev_or = [[0.895, 0.905, 0.62];       % vx=0.3
                           [0.39, 0.46, 0.46];      % vx=0.5
                           [0.245, 0.57, 0.61];  % vx=0.7
                           [0.54, 0.435, 0.5];  % vx=1
                           [0.455, 0.58 ,0.665];  % vx=-0.3
                           [0.03, 0.75, 0.115];  % vx=-0.5
                           [0.0, 0, 0.03];  % vx=-0.7
                           [0.05, 0.045, 0]   % vx=-1
                           ];                        % [m] Deviazione laterale ottenuta

Results.BGR4ET3.spost = [[0.7, 0.71, 0.79];          % vx=0.3
                          [1.20, 1.38, 1.405];         % vx=0.5
                          [1.61, 2.10, 1.88];     % vx=0.7
                          [1.5, 1.47, 1.603];     % vx=1
                          %[0., 0., 0.];     % vx=-0.3
                          [1.2086, 1.3350, 1.3561]
                          [1, 1.185, 1.285];     % vx=-0.5
                          [1.2, 1.23, 1.23];     % vx=-0.7
                          [1.25, 1.27, 1.3]      % vx=-1
                          ];                % [m] Spostamento longitudinale ottenuto

Results.BGR4ET3.dev_or = [[0.32, 0.37, 0.3];       % vx=0.3
                           [0.47, 0.44, 0.225];      % vx=0.5
                           [0.44, 0.2, 0.71];  % vx=0.7
                           [0.13, 0.195, 0.39];  % vx=1
                           %[0., 0., 0.];  % vx=-0.3
                           [0.2680, 0.2536, 0.3386]
                           [0.16, 0.24, 0.31];  % vx=-0.5
                           [0.216, 0.21, 0.31];  % vx=-0.7
                           [0.14, 0.12, 0.125]   % vx=-1
                           ];                        % [m] Deviazione laterale ottenuta

Results.EGREET3.spost = [%[0.29, 0.39, 0.42];          % vx=0.3
                          [1.3929, 1.4457, 1.4829]
                          [0.83, 0.885, 0.92];         % vx=0.5
                          [1.145, 1.24, 1.35];     % vx=0.7
                          [2.37, 2.35, 2.41];     % vx=1
                          [0.76, 0.8, 0.83];     % vx=-0.3
                          [1.235, 1.305, 1.36];     % vx=-0.5
                          [1.63, 1.695, 1.55];     % vx=-0.7
                          [1.78, 1.845, 1.96]      % vx=-1
                          ];                % [m] Spostamento longitudinale ottenuto

Results.EGREET3.dev_or = [[0., 0., 0.];       % vx=0.3
                          [0.215, 0.225, 0.21];      % vx=0.5
                          [0.59, 0.43, 0.35];  % vx=0.7
                          [0., 0.05, 0.23];  % vx=1
                          [0.06, 0.04, 0.02];  % vx=-0.3
                          [0.29, 0.26, 0.475];  % vx=-0.5
                          [0.545, 0.51, 0.32];  % vx=-0.7
                          [0.705, 0.635, 0.63]   % vx=-1
                           ];                        % [m] Deviazione laterale ottenuta

%% =================== CALCOLO DEI VARI TERMINI UTILI ================== %%
% Da qua in poi non serve toccare più nulla, i calcoli sono tutti automatizzati e
% si adattano da soli alle grandezze che inserisci.

%======== Riempimento della struttura e partizionamento dei dati =========%
Policy = struct();           % Inizializzazione delle strutture
for i =1:numel(Policy_name)
    Policy.(Policy_name{i}) = struct();                                    % Inizializzazione strutture policy
    Policy.(Policy_name{i}).vx_com = vx_com;                               % Salvo il comando di velocità
    Policy.(Policy_name{i}).spost_des = spost_des;                         % Salvo lo spostamento desiderato   
    Policy.(Policy_name{i}).durata_gradino = abs(Policy.(Policy_name{i}).spost_des ./ Policy.(Policy_name{i}).vx_com);  % Calcolo quanto deve durare il gradino da comandare
    vx_com_string = 'vx' + strrep(string(Policy.(Policy_name{i}).vx_com), '.', '_');   % Riporto gli elementi di vx come stringa per usarli come sottostrutture
    vx_com_string = strrep(vx_com_string, '-','_neg_');
    for j = 1:numel(vx_com_string)
        Policy.(Policy_name{i}).(vx_com_string(j)).Spost_act = struct();   % Inizializzo struttura in cui salvare i dati di spostamento longitudinale
        Policy.(Policy_name{i}).(vx_com_string(j)).Dev_or = struct();      % Inizializzo struttura in cui salvare i dati di deviazione laterale
        for k = 1:numel(index_test_results)                                % Riempio le strutture inizializzate
            Policy.(Policy_name{i}).(vx_com_string(j)).Spost_act.(index_test_results{k}) = Results.(Policy_name{i}).spost(j,k);
            Policy.(Policy_name{i}).(vx_com_string(j)).Dev_or.(index_test_results{k}) = Results.(Policy_name{i}).dev_or(j,k);
        end
    end
end

% ======================= Analisi statistica dati ======================= %
for i =1:numel(Policy_name)
    for j = 1:numel(vx_com_string)
    % Analisi dello spostamento longitudinale
    Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Spost.mean = mean(Results.(Policy_name{i}).spost(j,:));                                              % Calcolo media dei campioni
    Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Spost.err_mean = mean(Policy.(Policy_name{i}).spost_des - Results.(Policy_name{i}).spost(j,:));      % Calcolo la media dell'errore
    Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Spost.std_dev_err = std(Policy.(Policy_name{i}).spost_des - Results.(Policy_name{i}).spost(j,:));    % Calcolo la deviazione standard dell'errore
    % Analisi dello spostamento laterale
    Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Dev_or.mean = mean(Results.(Policy_name{i}).dev_or(j,:));              % Calcolo la media dell'errore
    Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Dev_or.std_dev_err = std(Results.(Policy_name{i}).dev_or(j,:));        % Calcolo la deviazione standard dell'errore
    end
end
%% ======================== PLOT DEI RISULTATI ========================= %%
perf_track = zeros(numel(Policy_name), numel(vx_com_string));
perf_dev_y = zeros(numel(Policy_name), numel(vx_com_string));
% ====================== Distribuzione dei campioni ===================== %
for i =1:numel(Policy_name)
    for j = 1:numel(vx_com_string)
        figure(i)
        sgtitle(Policy_name{i} + ": Distribuzione campioni spost x" )
        subplot(2,4,j)
        mean_res = mean(Results.(Policy_name{i}).spost(j,:));
        vx_mean = mean_res / Policy.(Policy_name{i}).durata_gradino(j);
        perf_track(i,j) = round(abs(vx_mean / vx_com(j) * 100), 3, "significant");
        bar(Results.(Policy_name{i}).spost(j,:), 'FaceColor',color(i)); hold on; yline(mean_res, '--r', 'LineWidth',2); yline(Policy.(Policy_name{i}).spost_des, 'g', 'LineWidth',2); hold off;
        title(strrep(vx_com_string(j), '_','.') + ' -> vx_{mean}= ' + vx_mean + '  \eta_x = ' + perf_track(i,j) + '%'); xlabel("# test"); ylabel("Spost. long [m]"); box on; grid on; 
        meanStr = sprintf('Media = %.2f [m]', mean_res); legend(["Dati", meanStr, "Valore Desiderato"])
    end
end

for i =1:numel(Policy_name)
    for j = 1:numel(vx_com_string)
        figure(i+numel(Policy_name))
        sgtitle(Policy_name{i} + ": Distribuzione campioni deviazione y")
        subplot(2,4,j)
        mean_res = mean(Results.(Policy_name{i}).dev_or(j,:));
        perf_dev_y(i,j) = round((- abs(mean_res / spost_des) * 100), 3, "significant");
        bar(Results.(Policy_name{i}).dev_or(j,:), 'FaceColor',color(i)); hold on; yline(mean_res, '--r', 'LineWidth',2);hold off;
        title (strrep(vx_com_string(j), '_','.')+ '  \eta_y = ' + perf_dev_y(i,j) + '%'); xlabel("# test"); ylabel("Spost. long [m]"); box on; grid on; 
        meanStr = sprintf('Media = %.2f [m]', mean_res); legend(["Dati", meanStr])

    end
end

%  ============================== Errore medio ========================== %
for i =1:numel(Policy_name)
    figure(2*numel(Policy_name)+1);
    sgtitle("Errore medio spost x [m]")
    subplot(1,numel(Policy_name),i)
    Mean = zeros(1,8);
    for j = 1:numel(vx_com_string)
        Mean(j) = Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Spost.err_mean;
    end
    bar(Mean, 'FaceColor',color(i)); hold on; yline(mean(Mean, 'all'), '--r', 'LineWidth',2); hold off
    xticklabels(strrep(vx_com_string, '_','.')); ylabel("Policy: " + (Policy_name{i})); 
    meantotStr = sprintf('Media totale = %.2f [m]',mean(Mean)); legend(["Errore medio" , meantotStr]); box on; grid on;
end

for i =1:numel(Policy_name)
    figure(2*numel(Policy_name)+2);
    sgtitle("Media della deviazione y [m]")
    subplot(1,numel(Policy_name),i)
    Mean = zeros(1,8);
    for j = 1:numel(vx_com_string)
        Mean(j) = Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Dev_or.mean;
    end
    bar(Mean, 'FaceColor',color(i)); hold on; yline(mean(Mean, 'all'), '--r', 'LineWidth',2); hold off
    xticklabels(strrep(vx_com_string, '_','.')); ylabel("Policy: " + (Policy_name{i})); 
    meantotStr = sprintf('Media totale = %.2f [m]',mean(Mean)); legend(["Errore medio" , meantotStr]); box on; grid on;
end

%  ========================== Deviazioni Standard ======================= %
for i =1:numel(Policy_name)
    figure(2*numel(Policy_name)+3);
    sgtitle("Dev Std errore di spost x [m]")
    subplot(1,numel(Policy_name),i)
    Stdev = zeros(1,8);
    for j = 1:numel(vx_com_string)
        Stdev(j) = Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Spost.std_dev_err;
    end
    bar(Stdev, 'FaceColor',color(i)); hold on; yline(mean(Stdev, 'all'), '--r', 'LineWidth',2); hold off
    xticklabels(strrep(vx_com_string, '_','.')); ylabel("Policy: " + (Policy_name{i})); 
    meanStdStr = sprintf('Std. Dev. media = %.2f [m]',mean(Stdev)); legend(["Std. Dev." , meanStdStr]); box on; grid on;
end

for i =1:numel(Policy_name)
    figure(2*numel(Policy_name)+4);
    sgtitle("Dev Std della deviazione y [m]")
    subplot(1,numel(Policy_name),i)
    Stdev = zeros(1,8);
    for j = 1:numel(vx_com_string)
        Stdev(j) = Policy.(Policy_name{i}).(vx_com_string(j)).Analisi.Dev_or.std_dev_err;
    end
    bar(Stdev, 'FaceColor',color(i)); hold on; yline(mean(Stdev, 'all'), '--r', 'LineWidth',2); hold off
    xticklabels(strrep(vx_com_string, '_','.')); ylabel("Policy: " + (Policy_name{i})); 
    meanStdStr = sprintf('Std. Dev. media = %.2f [m]',mean(Stdev)); legend(["Std. Dev." , meanStdStr]); box on; grid on;
end

%  ============================== Performance =========================== %
mean_perf_track = mean(perf_track, 2);
mean_perf_dev_y = mean(perf_dev_y, 2);
perf_tot = perf_track + perf_dev_y;
mean_perf_tot = mean(perf_tot, 2);

figure(2*numel(Policy_name)+5);
sgtitle("Performance complessive");
subplot(1,3,1)
    b1 = bar(mean_perf_track);
    xticklabels(Policy_name); box on; grid on; ylabel("[%]"); title("\eta_x")
    % Impostazione dei colori delle barre
    for k = 1:numel(Policy_name)
        b1.FaceColor = 'flat';
        b1.CData(k,:) = sscanf(color(k).char, '#%2x%2x%2x', [1 3])/255;  % converte i colori da HEX a RGB
    end
subplot(1,3,2)
    b2 = bar(mean_perf_dev_y);
    xticklabels(Policy_name); box on; grid on; ylabel("[%]"); title("\eta_y")
    % Impostazione dei colori delle barre
    for k = 1:numel(Policy_name)
        b2.FaceColor = 'flat';
        b2.CData(k,:) = sscanf(color(k).char, '#%2x%2x%2x', [1 3])/255;  % converte i colori da HEX a RGB
    end
subplot(1,3,3)
    b3 = bar(mean_perf_tot);
    xticklabels(Policy_name); box on; grid on; ylabel("[%]"); title("\eta_{tot}")
    % Impostazione dei colori delle barre
    for k = 1:numel(Policy_name)
        b3.FaceColor = 'flat';
        b3.CData(k,:) = sscanf(color(k).char, '#%2x%2x%2x', [1 3])/255;  % converte i colori da HEX a RGB
    end
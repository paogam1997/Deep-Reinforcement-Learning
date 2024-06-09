function [complete_delta, delta_, repetition] = FinalMaxDelta(parametro, episode_lenght)
if parametro.frequency_interval == 0
    parametro.frequency_interval = 999999;  % Riconoscibile in modo da non fare il conto se on_interval non c'è
end
int = floor(episode_lenght /parametro.frequency_interval);
delta = zeros(int+1, 1);
% Primo indice corrisponde a condizione on_reset che sta nella prima riga
% di .dr_info
if parametro.dr_info(1, 2) == "additive" && parametro.dr_info(1, 3) == "gaussian"
    if parametro.type == "angle"
        delta(1) = VarToDelta_angoli(parametro.distribution_parameters(1,1));
    elseif parametro.type == "general"
        delta(1) = VarToDelta_generale(parametro.distribution_parameters(1,1));
    end
elseif parametro.dr_info(1, 2) == "scaling" && parametro.dr_info(1, 3) == "uniform"
    delta(1) = parametro.distribution_parameters(1,2); % MAX
elseif parametro.dr_info(1, 2) == "additive" && parametro.dr_info(1, 3) == "uniform"
    if parametro.distribution_parameters(1,1) == - parametro.distribution_parameters(1,2)  % usi un intervallo simmetrico
        delta(1) = parametro.distribution_parameters(1,2);
    else
        if abs(parametro.distribution_parameters(1,1)) > abs(parametro.distribution_parameters(1,2))
            segno = sign(parametro.distribution_parameters(1,1));
        else
            segno = sign(parametro.distribution_parameters(1,2));
        end
        delta(1) = segno* max(abs(parametro.distribution_parameters(1,1)),abs(parametro.distribution_parameters(1,2)));
    end
elseif parametro.dr_info(1, 2) == "direct"
    delta(1) = 0;
end


% Indici successivi corrispondono a condizione on_interval che sta nella seconda riga
% di .dr_info
for j = 1:int
    if parametro.dr_info(2, 2) == "additive" && parametro.dr_info(2, 3) == "gaussian"
        if parametro.type == "angle"
            delta(j+1) = VarToDelta_angoli(parametro.distribution_parameters(2,1));
        elseif parametro.type == "general"
            delta(j+1) = VarToDelta_generale(parametro.distribution_parameters(2,1));
        end
    elseif parametro.dr_info(2, 2) == "scaling" && parametro.dr_info(2, 3) == "uniform"
        delta(j+1) = parametro.distribution_parameters(2,2);
    elseif parametro.dr_info(2, 2) == "additive" && parametro.dr_info(2, 3) == "uniform"
        if parametro.distribution_parameters(2,1) == - parametro.distribution_parameters(2,2)  % usi un intervallo simmetrico
            delta(j+1) = parametro.distribution_parameters(2,2);
        else
            if abs(parametro.distribution_parameters(2,1)) > abs(parametro.distribution_parameters(2,2))
                segno = sign(parametro.distribution_parameters(2,1));
            else
                segno = sign(parametro.distribution_parameters(2,2));
            end
            delta(j+1) = segno* max(abs(parametro.distribution_parameters(2,1)),abs(parametro.distribution_parameters(2,2)));
        end
    elseif parametro.dr_info(2, 2) == "direct"
        delta(j+1) = 0;
    end
end

if parametro.dr_info(2, 2) == "None" % Non c'è una randomizzazione on_interval
    complete_delta = delta(1);
    int = 0;
else
    if parametro.dr_info(1, 2) == "additive" && (parametro.dr_info(1, 3) == "gaussian" || parametro.dr_info(1, 3) == "uniform")

        if parametro.dr_info(2, 2) == "additive" && (parametro.dr_info(2, 3) == "gaussian" || parametro.dr_info(2, 3) == "uniform")
            complete_delta = sum(delta, "all");
        elseif parametro.dr_info(2, 2) == "scaling" && parametro.dr_info(2, 3) == "uniform"
            complete_delta = delta(1)*delta(2:end); % sommi la perturbazione di on_reset e quelle on_interval scalano sul pertubato
        end
    elseif parametro.dr_info(1, 2) == "scaling" && parametro.dr_info(1, 3) == "uniform"
        if parametro.dr_info(2, 2) == "scaling" && parametro.dr_info(2, 3) == "uniform"
            complete_delta = prod(delta);
        elseif parametro.dr_info(2, 2) == "additive" && (parametro.dr_info(2, 3) == "gaussian" || parametro.dr_info(2, 3) == "uniform")
            complete_delta = delta(1)+delta(2:end);
        end
    elseif parametro.dr_info(1, 2) == "direct"
        if parametro.dr_info(2, 2) == "scaling" && parametro.dr_info(2, 3) == "uniform"
            complete_delta = prod(delta(2:end));
        elseif parametro.dr_info(2, 2) == "additive" && (parametro.dr_info(2, 3) == "gaussian" || parametro.dr_info(2, 3) == "uniform")
            complete_delta = sum(delta(2:end));
        end
        % elseif parametro.dr_info(1, 2) == "additive" && parametro.dr_info(1, 3) == "uniform"
        %     complete_delta = sum(delta, "all");
    end

end
delta_ = delta(:,:);
repetition = int;
end
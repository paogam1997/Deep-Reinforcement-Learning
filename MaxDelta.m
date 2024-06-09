function max_variation = MaxDelta(parametro, method)
    if method == "on_reset"
        i = 1;
    elseif method == "on_interval"
        i = 2;
    end

    if parametro.dr_info(1, 2) == "additive" && parametro.dr_info(1, 3) == "gaussian"
        if parametro.type == "angle"
            max_variation = VarToDelta_angoli(parametro.distribution_parameters(i,1));
        elseif parametro.type == "general"
            max_variation = VarToDelta_generale(parametro.distribution_parameters(i,1));
        end
    elseif parametro.dr_info(1, 2) == "scaling" && parametro.dr_info(1, 3) == "uniform"
            max_variation = parametro.distribution_parameters(i,2);
    elseif parametro.dr_info(1, 2) == "additive" && parametro.dr_info(1, 3) == "uniform"
        if parametro.distribution_parameters(i,1) == - parametro.distribution_parameters(i,2)  % usi un intervallo simmetrico
            max_variation = parametro.distribution_parameters(i,2);
        else
            if abs(parametro.distribution_parameters(i,1)) > abs(parametro.distribution_parameters(i,2))
                segno = sign(parametro.distribution_parameters(i,1));
            else
                segno = sign(parametro.distribution_parameters(i,2));
            end
            max_variation = segno* max(abs(parametro.distribution_parameters(i,1)),abs(parametro.distribution_parameters(i,2)));
        end
        elseif parametro.dr_info(1, 2) == "direct"
            max_variation = parametro.distribution_parameters(i,1);
    end
end
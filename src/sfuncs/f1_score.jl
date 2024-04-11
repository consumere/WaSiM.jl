# 
function f1_score(data::DataFrame, 
        sim_col::Union{Symbol,Int64} = 1, 
        obs_col::Union{Symbol,Int64} = 2; threshold = 0.5, minobs=true)
        # Extracting simulation and observation columns
        #sim_values = data[:, sim_col]
        #obs_values = data[:, obs_col]
        sim_nm=only(names(select(data, sim_col)))
        obs_nm=only(names(select(data, obs_col)))
        @info "simcol is $sim_nm, obscol is $obs_nm"
        sim_values = select(data, sim_col)|>Matrix|>vec
        obs_values = select(data, obs_col)|>Matrix|>vec

        if minobs
            threshold = minimum(obs_values)
            @info "threshold is set to minimum of observations: $threshold"
        else
            # Threshold for classification (you may need to adjust this based on your problem)
            threshold = threshold
        end

        # Convert to binary classification (1 if greater than threshold, 0 otherwise)
        sim_binary = sim_values .> threshold
        obs_binary = obs_values .> threshold

        # True Positives, False Positives, True Negatives, False Negatives
        tp = sum(sim_binary .& obs_binary)
        fp = sum(sim_binary .& .!obs_binary)
        tn = sum(.!sim_binary .& .!obs_binary)
        fn = sum(.!sim_binary .& obs_binary)

        # Precision, Recall, and F1 Score
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        
        # F1 Score is the harmonic mean of precision and recall
        f1_score = 2 * (precision * recall) / (precision + recall)

        return f1_score
    end

    """
    recursive rmeq; keeps otherdata
    """
    
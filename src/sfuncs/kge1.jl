# 
function kge1(simulations, evaluation)
        r = cor(simulations, evaluation)
        α = std(simulations) / std(evaluation)
        β = mean(simulations) / mean(evaluation)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    
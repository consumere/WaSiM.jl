# 
function nse(df::DataFrame)
        simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
        return (1 - (sum((simulated .- observed).^2) / sum((observed .- mean(observed)).^2)))
    end

    
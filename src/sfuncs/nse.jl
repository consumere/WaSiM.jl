# 
function nse(df::DataFrame)
        simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
        #observed, simulated = df[:,6],df[:,5]
        return (1 - (sum((simulated .- observed).^2) / sum((observed .- mean(observed)).^2)))
    end

    """
    kge as in Gupta et al., 2009
    https://doi.org/10.1016/j.jhydrol.2009.08.003
    """
    
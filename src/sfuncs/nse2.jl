# 
function nse2(df::DataFrame)
        observed, simulated = df[:,6],df[:,5]
        return (1 - (sum((simulated .- observed).^2) / sum((observed .- mean(observed)).^2)))
    end
    
    """
    reads recursively
    """
    
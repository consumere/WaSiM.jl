# 
function vef2(df::DataFrame)
        obs, sim = df[:,6],df[:,5]
        return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
    end

    """
    reads recursively from qout and calulates KGE and NSE
    sorted by NSE
    """
    
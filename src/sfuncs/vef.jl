# 
function vef(x::String)
        ms=["-9999","lin","log"]
        df::DataFrame = CSV.read(x,DataFrame,
        missingstring=ms,
        types = Float64,
        delim="\t",
        silencewarnings=true,
        normalizenames=true) |> dropmissing
        #drop=(i, nm) -> i == 4) |> dropmissing
        
        obs, sim = df[:,6],df[:,5]
        
        return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
    end

    """
    selects columns from vector of dfs ...
    """
    
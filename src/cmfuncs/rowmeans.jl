# 
function rowmeans(df::DataFrame)
        return DataFrames.combine(df, names(df) .=> mean)
    end
    
    
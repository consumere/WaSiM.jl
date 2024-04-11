# 
function rowsums(df::DataFrame)
        return DataFrames.combine(df, names(df) .=> sum)
    end
    
    
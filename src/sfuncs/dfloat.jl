
"""
tries to parse all columns to float, except of date.
"""
function dfloat(x::DataFrame)
    df = copy(x)
    for colname in names(df)
        if eltype(df[!, colname]) != Date && eltype(df[!, colname]) == String
            #&& eltype(df[!, colname]) != Float64
            df[!, colname] .= tryparse.(Float64, df[!, colname])
        end
    end
    return df
end    

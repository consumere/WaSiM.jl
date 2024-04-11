# 
function lplotf(df::String)
        df=readdf(df)
        nm = propertynames(df)[1:end-1];
        o = collect(DataFrames.metadata(df))[1][2] |>basename
        @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=o)     
    end

    
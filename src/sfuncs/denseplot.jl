# 
function denseplot(df::String)
        df=readdf(df)
        s = propertynames(df)[Not(end)]
        @df df density(cols(s), legend = :topright)
    end


    
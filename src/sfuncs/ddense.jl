# 
function ddense(df::DataFrame)
        s = Symbol.(filter(x->!occursin(r"date|year"i,x),names(df)))
        @df df density(cols(s), legend = :topright)
    end

    
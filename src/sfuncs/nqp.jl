# 
function nqp(a::Regex,b::Regex;)
        a = waread(a)
        b = waread(b)
        col = ncol(a)-1
    
        a = a[!,Cols(col,:date)]
        b = b[!,Cols(col,:date)]  
    
        df = mall(a,b)
        #qplot(x::Vector{Float64},y::Vector{Float64})
        return qplot(df)
    end

    
# 
function kge(df::DataFrame;verbose=false)
        if (any(x->occursin("year|date|month",x),names(df)))
            ln = Symbol.(filter(x->!occursin(r"date|year|month"i,x),names(df)))
            df = select(df, ln)
        end
        simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
        r = cor(observed, simulated)
        α = std(simulated) / std(observed)
        β = mean(simulated) / mean(observed)
        kge = 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
        if verbose
            println("simulated: ",names(df)[1])
            println("observed: ",names(df)[2])
            #(r, α, β) as per `Gupta et al., 2009
            println("returning r, α, β, kge")
            return DataFrame(r=r, α=α, β=β, kge=kge)
        else
            return kge
        end
        
    end

    
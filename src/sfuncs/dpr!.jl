# 
function dpr!(x::DataFrame)
        """
        correlation plots on dataframe
        """
        df = copy(x)
        if any(map(x->occursin("year",x),names(df)))
            df = df[!,Not(:year)]
        end
        
        if any(map(x->occursin("month",x),names(df)))
            df = df[!,Not(:month)]
        end
        
        if propertynames(df)[end]!=:date
            df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
        end
        
        dropmissing!(df)
        
        v = (names(df[!,1:2]))
        a = reshape(v, 1, 2)|>reverse #sim,obs
        Plots.plot(df.date,[df[!,2], df[!,1]], 
        label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
        r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
        # annotate!(last(df.date), 0.85*maximum(df[!,1]),
        # text("R² = $r2", 10, :black, :right))
        dropmissing!(df)
        kge = round(kge2(df[!,2], df[!,1]), digits=2)
        nse_value = round(nse(df[!,1], df[!,2]), digits=2)
        ve = round(vef(df[!,2], df[!,1]), digits=2)
        annotate!(
            last(df.date), 0.95*maximum(df[!,1]),
        text("KGE = $kge\nNSE = $nse_value\nR² = $r2\nVE = $ve\n", 
        10, :black, :right)
        )
    end
    
    
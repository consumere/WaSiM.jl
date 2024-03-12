# 
function dprbig(x::Regex)
        df = globdf(x)|>first|>waread
        df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
        v = (names(df[!,1:2]))
        a = reshape(v, 1, 2)|>reverse #sim,obs
        Plots.plot(df.date,[df[!,2], df[!,1]], 
        xlabel="Date", ylabel="[mm/day]",
        legend = :topleft,
        size=(1200,800)
        )
        r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
        kge = round(kge2(df[!,2], df[!,1]), digits=2)
        nse_value = round(nse(df[!,1], df[!,2]), digits=2)
        annotate!(
            last(df.date), 0.95*maximum(df[!,1]),
        text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
        )
    end

    #
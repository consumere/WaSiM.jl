# 
function dpr(a::Regex,b::Regex)
        """
        correlation plots on dataframe
        """
        a = waread(a)
        b = waread(b)
        # colA = ncol(a)-1
        # colB = ncol(b)-1

        # a = a[!,Cols(colA,:date)]
        # b = b[!,Cols(colB,:date)]
        
        a = a[!,Cols(1,:date)]
        b = b[!,Cols(1,:date)] 

        df = mall(a,b)
        dropmissing!(df)
        
        df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
        v = (names(df[!,1:2]))
        a = reshape(v, 1, 2)

        Plots.plot(df.date,[df[!,1], df[!,2]],  label=a, 
        #    seriestype = :bar,
            xlabel="Date", ylabel="[mm/day]",
            legend = :topleft)

        r2 = round(cor(df[!,1], df[!,2])^2, digits=2)
        kge = round(kge2(df[!,2], df[!,1]), digits=2)
        nse_value = round(nse(df[!,1], df[!,2]), digits=2)
        ve = round(vef(df[!,2], df[!,1]), digits=2)

        annotate!(
            last(df.date), 0.95*maximum(df[!,1]),
        text("KGE = $kge\nNSE = $nse_value\nR² = $r2\nVE = $ve\n", 
        10, :black, :right)
        )
    end

    
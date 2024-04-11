# 
function monp(lk::DataFrame;doleg=true)
        
        df = monmean(lk)
        #str = [ @sprintf("%02i", x) for x in (df.month) ];
        #mn = [ monthname(x) for x in (df.month) ]
        mn = [ monthabbr(x) for x in (df.month) ]
        ln = Symbol.(filter(x->!occursin(r"date|month|year",x),names(df)))
        Plots.plot()
    
        if doleg
            for (i, col) in enumerate(ln)
                Plots.plot!(mn, df[!,col], fillalpha=0.75, 
                #linestyle=linestyles[i], 
                linestyle=:auto, 
                linewidth=1.25, 
                #label=[string(col)],
                label=string(col),
                legend=true)
            end
    
            return Plots.plot!() 
        end 
    
        for (i, col) in enumerate(ln)
            Plots.plot!(mn, df[!,col], fillalpha=0.75, 
            linestyle=:auto, 
            linewidth=1.25, 
            legend=false)
        end
        return Plots.plot!() 
    end

    
# 
function qplot(df1::DataFrame,df2::DataFrame;col1=1,col2=1)
        if any(map(x->contains(x,"date"),names(df1)))
            df1 = df1[!,Not(:date)]
        end

        if any(map(x->contains(x,"date"),names(df2)))
            df2 = df2[!,Not(:date)]
        end

        try
            df1 = select!(df1, col1)
            df2 = select!(df2, col2)
        catch
            @error "col not found!"
            return
        end       

        # if ncol(df)>2
        #     df = df[!,1:2]
        # end

        r2 = round(cor(df1[!,1], df2[!,1])^2, digits=3)
        
        p = qqplot(df1[!,1], df2[!,1], 
            #title = "R² = "*string(ti),
            qqline = :fit)
            #color = :grays) # erstellt ein QQ-Diagramm <- black
        xlabel!(p,names(df1)[1])
        ylabel!(p,names(df2)[1])
        annotate!(p,:bottomright, text("R² = "*string(r2), :black))
                
    end

    """
    greps first from current dir Regex and copies to clipboard
    """
    
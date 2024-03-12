# 
function cmplot(;temp::Union{Regex,String}=r"^so_temper",prec::Union{Regex,String}=r"^so_prec",col=1)
        #temp,prec,col = (r"^so_temper",r"^so_pre",1)
        #eig. braucht man das col argument bei so_ input nicht.

        prec = waread2(prec;silencewarnings=true,limit=10^6)
        
        yrs = year.(prec.date)|>unique|>length
        prec = monsum(prec)
        if ncol(prec) > 2
            select!(prec, Not(:month))
            prec = prec[:,col]
            precvec = vec(Matrix(prec))
        else
            precvec = vec(Matrix(select(prec, Not(:month))))
        end
                
        precvec = precvec ./ yrs
        
        temp = waread2(temp;
            silencewarnings=true,    
            #maxwarnings=1,
            limit=10^6)|>monmean
        if ncol(temp) > 2
            select!(temp, Not(:month))
            temp = temp[:,col]
            tempvec = vec(Matrix(temp))
        else
            tempvec = vec(Matrix(select(temp, Not(:month))))
        end
        
        month_abbr = ["Jan", "Feb", "Mär", "Apr", 
            "Mai", "Jun", "Jul", "Aug", "Sep", 
                "Okt", "Nov", "Dez"];
        
        #Plots.showtheme(:dao)
        
               
         #;family="Computer Modern"
        #Plots.theme(:ggplot2) #passt.
        #Plots.theme(:dao) #gibt mir falsche twinx margin.
        #Plots.theme(:sheet)
        #Plots.theme(:dracula)
        Plots.theme(:wong2)
        #Plots.showtheme(:dracula)
        #Plots.showtheme(:dracula)
        
        p1 = Plots.bar(prec.month, precvec, 
            color=:cornflowerblue, 
            guidefontfamily="Computer Modern",
            tickfontfamily="Computer Modern",
            #xlabel="", 
            ylims = (0.0, maximum(precvec) + 5.0),
            xflip=false,
            ylabel="Niederschlag [mm]", 
            legend=false, yflip=true,
            tick_direction = :out
            );
        xticks!(1:12, month_abbr)

        for i in prec.month
            val = round(precvec[i]; digits=1)
            annotate!(i, precvec[i], 
            Plots.text("$(val)",8,:bottom; family="Computer Modern"))
        end
        
        #font(family="serif", halign=:center, rotation=45.0)
        # val = round(maximum(tempvec); sigdigits=2)
        # i = findmax(tempvec)[2]
        # Plots.annotate!(i, tempvec[i], 
        #             text("max: $(val) °C",10,:top))
        
        #p2 = twinx(); #s.u.
        
        ann2 = map(x->
            Plots.text(
            string.(round(x; digits=1))*"°", 
                8,
                :left, 
                #:red),
                :black;
                family="Computer Modern"),
                #family="sans-serif"), #maby
                tempvec)
        
        plot!(twinx(), tempvec, #xlabel="", 
            ylabel="Temperatur [°C]", 
            guidefontfamily="Computer Modern",
            tickfontfamily="Computer Modern",
            color=:coral2,
            ylims = (minimum(tempvec) - 1.0 , 
                maximum(tempvec) + 1.0),
            annotations = (temp.month .+ 0.125, 
                tempvec .+ 0.1, ann2, :center),
            label=false, 
            linestyle = :dashdot,
            linewidth = 2,
            right_margin = 5mm,
            tick_direction = :out            
            );      

        #@show plotattr(:Axis)

        return p1
    end

    """
    selects from dataframe date an column
    selt(x::DataFrame, col::Any;dtcol=:date)
    """
    
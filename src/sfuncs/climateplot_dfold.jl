# 
function climateplot_dfold(temp::DataFrame,prec::DataFrame;col::AbstractString)
        #col = propertynames(temp)[col]
        #temp = t
        
        #col = first(Symbol.(filter(x->occursin(r"$col"i,x),names(temp))))
        #prec = pr
        yrs = year.(prec.date)|>unique|>length
        prec = monsum(prec)
        precvec = vec(Matrix(select(prec, col)))
        #precvec = vec(Matrix(select(prec, Not(:month))))
        precvec = precvec ./ yrs
        
        temp = (temp)|>monmean
        tempvec = vec(Matrix(select(temp, col)))
        #tempvec = vec(Matrix(select(temp, Not(:month))))
        month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
        p1 = Plots.plot(temp.month,tempvec, xlabel="", 
            ylabel="Temperature [°C]", color=:coral2,
            yflip = false,
            label=false, linewidth=3)
        # # Add annotations for temperature values
        for i in 1:length(tempvec)
            val = round(tempvec[i]; sigdigits=1)
            annotate!(i, tempvec[i], text("$(val)",7, :center))
        end
        
        #twinx()
        #prec.month, 
        Plots.bar!(precvec, color=:cornflowerblue, 
        xlabel="",     #xlabel="Months", 
        xflip=false,
        ylabel="Precipitation [mm]", 
        legend=false    #, yflip=true
        )
        xticks!(1:12, month_abbr)
        for i in prec.month
            val = round(precvec[i]; digits=1)
            annotate!(i, precvec[i], text("$(val)",8, :bottom))
        end
        return p1
    end

    
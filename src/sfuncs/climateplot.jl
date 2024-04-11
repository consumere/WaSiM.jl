# 
function climateplot(temp::Union{DataFrame,String},prec::Union{DataFrame,String},col::Union{Symbol,Int,String} = 1)
        if typeof(temp)==String
            temp = waread(temp)
        end
        if typeof(prec)==String
            prec = waread(prec)
        end
        yrs = year.(prec.date)|>unique|>length
        prec = monsum(prec)
        precvec = vec(Matrix(select(prec, col)))
        #precvec = vec(Matrix(select(prec, Not(:month))))
        precvec = precvec ./ yrs
        
        temp = (temp)|>monmean
        tempvec = vec(Matrix(select(temp, col)))
        month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
        p1 = Plots.bar(prec.month, precvec, color=:cornflowerblue, 
            xlabel="", 
            #xlabel="Months", 
            xflip=false,
            ylabel="Precipitation [mm]", 
            #ylims = (0, Int64(ceil(maximum(precvec) + 20))),
            ylims = (0, Int64(ceil(maximum(precvec)*1.1))),
            legend=false, yflip=true);

        #ylims!(0, Int64(ceil(maximum(tempvec) + 20)))
        
        xticks!(1:12, month_abbr)
        for i in prec.month
            val = round(precvec[i]; digits=1)
            annotate!(i, precvec[i], text("$(val)",8, :bottom))
        end
        p2 = twinx();
        ann2 = map(x->
        text(
        string.(round(x; digits=1)) , 7 ,:center, :red
        )        ,tempvec)
        plot!(p2, tempvec, xlabel="", 
            ylabel="Temperature [°C]", color=:coral2,
            annotations = (temp.month,tempvec, ann2),
            #fontsize=7,
            label=false, linewidth=3);
            
        # # # Add annotations for temperature values
        # for i in 1:length(tempvec)
        #     val = round(tempvec[i]; sigdigits=3)
        #     annotate!(i, tempvec[i], text("$(val)",7, :center))
        # end
        return p1
    end

    """
    ws. prc and temp tauschen un opacitiy einstellen....
    ,col::String name of 
    twinx() dreht komplett alles.
    """
    
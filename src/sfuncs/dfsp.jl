# 
function dfsp(dfs::Vector{DataFrame};save="")
        "plots and adds"
        df = dfs[1]
        s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
        p = @df df Plots.plot(:date,
                cols(s),
                legend = false)    #legend = :bottom)
        for i in 2:length(dfs)
            nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
            println("adding $nm")
            s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
            @df dfs[i] Plots.plot!(:date,cols(s),
            label="$nm") #geht, wenn oben legend true ist.
        end
        return p
        if !isempty(save) 
            Plots.savefig(p,save*".png")
            printstyled("$save saved as $save*.png! \n",color=:green)
        end
    end

    
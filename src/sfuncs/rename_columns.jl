# 
function rename_columns(x::DataFrame, name_mapping::DataFrame)
        df = copy(x)
        v = map(x->Grep.grep(x,names(df)),Regex.(string.(name_mapping.sim)))
        for (i, name) in enumerate(reduce(vcat, v))
            if occursin(string.(name_mapping.sim[i]),name)
                println("renaming ",name, " <-> ",name_mapping.sim[i], " => " ,name_mapping.name[i])
                rename!(df, name => name_mapping.name[i])
            end
        end
        return(df)
    end


    """
    adds trendlines to plot
    tline(df::DataFrame, date_col::Symbol;lab=true)
    lab sets the trendline equation as label
    """
    
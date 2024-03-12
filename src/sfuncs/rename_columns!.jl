# 
function rename_columns!(df::DataFrame, name_mapping::DataFrame)
        v = map(x->Grep.grep(x,names(df)),Regex.(string.(name_mapping.sim)))
        for (i, name) in enumerate(reduce(vcat, v))
            if occursin(string.(name_mapping.sim[i]),name)
                println("renaming ",name, " <-> ",name_mapping.sim[i], " => " ,name_mapping.name[i])
                rename!(df, name => name_mapping.name[i])
            end
        end
    end

    """
    rename_columns(df::DataFrame, name_mapping::DataFrame)
    """
    
# 
function dfbar(df::DataFrame)
        ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        if (any(x->occursin("year",x),names(df)))
            s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
            #@df df Plots.bar(:year,cols(s),legend = :topright, title=ti)
            @df df groupedbar(df.year,cols(s), legend = :outertopright, title=ti)
        else    
        s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
        #@df df Plots.bar(:date,cols(s),legend = :topright, title=ti)
        @df df groupedbar(df.date,cols(s), legend = :outertopright, title=ti)
        end
    end

    """
    return 1st DFs of qgko see also qbb for rglob
    """
    
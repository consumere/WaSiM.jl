# 
function dfl(x::DataFrame)
        "selects first match and plots in log y-axis..."
        df=x
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw"" 
        end
        if (any(x->occursin("year",x),names(df)))
            s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
            @df df Plots.plot(:year,cols(s),yaxis=:log,legend = :topright, title=ti)
        else   
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot(:date,cols(s),yaxis=:log10, 
            legend = :topright, title=ti)
            end
    end
    
    """
    rglob qgko
    return a vector of DFs
    """
    
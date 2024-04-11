# 
function dfl!(regex::Union{Regex,String})
        df=readf(regex)
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
        ti = raw"" 
        end
        println("adding $ti")
        if (any(x->occursin("year",x),names(df)))
            s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
            @df df Plots.plot!(:year,cols(s),yaxis=:log,legend = :topright)
            Plots.annotate!([(20,5,text(ti, 12, :left, :top, :green))])
        else   
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot!(:date,cols(s),yaxis=:log, legend = :topright)
            Plots.annotate!([(20,5,text(ti, 12, :left, :top, :green))])
            end
    end

    
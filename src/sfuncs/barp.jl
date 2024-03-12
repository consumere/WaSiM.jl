# 
function barp(x::DataFrame)
        "with DataFrame input"
            df = x
                ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
            if any(x->occursin("year",x),names(df))
                ln = Symbol.(filter(x->!occursin("year",x),names(df)))
                @df df Plots.plot(:year,
                    cols(ln),
                    legend = :topright, 
                    title=ti,
                    seriestype=:bar) #color=:lightrainbow
            elseif any(x->occursin("month",x),names(df))
                ln = Symbol.(filter(x->!occursin("month",x),names(df)))
                @df df Plots.plot(:month,
                    cols(ln),
                    legend = :topright, 
                    title=ti,
                    seriestype=:bar)
            elseif (
                any(x->occursin("month",x),names(df)) & 
                any(x->occursin("year",x),names(df))            
                )
                ln = (filter(x->!occursin("month",x),names(df)))
                ln = Symbol.(filter(x->!occursin("year",x),ln))
                @df df Plots.plot(:month,
                    cols(ln),
                    legend = :topright, 
                    title=ti,
                    seriestype=:bar)
            else
                dfp(df)        
            end
    end

    
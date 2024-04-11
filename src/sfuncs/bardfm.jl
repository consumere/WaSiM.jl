# 
function bardfm(x::DataFrame)
        "with DataFrame input"
            df = x
            y = filter(x->!occursin("date",x),names(df))
            s = map(y -> Symbol(y),y)
                ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
            df[!, :year] = year.(df[!,:date]);
            df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
            @df df_yearsum Plots.plot(:year,
                cols(s),
                legend = :topright, 
                title=ti,
                seriestype=:bar)
    end

    
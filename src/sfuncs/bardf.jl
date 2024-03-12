# 
function bardf(x::DataFrame;leg=:topright)
            df = copy(x)
            y = filter(x->!occursin(r"date|year|month",x),names(df))
            s = map(y -> Symbol(y),y)
                ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
            #df[!, :year] = year.(df[!,:date]);
            #df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        nm = filter(x->occursin(r"date|year|month",x),
            names(df))|>first
        
        @df df Plots.plot(select(df,Symbol(nm))|>Matrix,
                cols(s),
                legend = leg, 
                title = ti,
                seriestype=:bar)
    end

    
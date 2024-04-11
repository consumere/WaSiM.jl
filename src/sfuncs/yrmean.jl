# 
function yrmean(x::DataFrame)
        df = copy(x)
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
        return(df_yearsum)
    end

    

    
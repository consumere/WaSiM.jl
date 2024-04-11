# 
function monsum(x::DataFrame)
        df = copy(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(groupby(df, :month), y .=> sum .=> y);
        return(df_monthsum)
    end

    
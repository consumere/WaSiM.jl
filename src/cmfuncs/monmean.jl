# 
function monmean(x::DataFrame)
        df = x
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(groupby(df, :month), y .=> mean .=> y);
        return(df_monthsum)
    end


    
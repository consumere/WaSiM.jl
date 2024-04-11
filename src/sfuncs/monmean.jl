# 
function monmean(x::DataFrame)
        df = copy(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        dmean = DataFrames.combine(groupby(df, :month), y .=> mean .=> y);
        return(dmean)
    end

    """
    running monthly mean...
    rmm(x::DataFrame;fun=mean)
    can also be sum. 
    dfr(r"qbas")|>i->rmm(i;fun=sum)|>cmk.dfp
    """
    
# 
function qrtr(pt::Union{String,DataFrame};fun=sum,agg=quarterofyear)

        if pt isa String
            x = waread2(pt)
        else
            x = pt
        end
        
        df = copy(x)
        y = filter(x->!occursin("date",x), names(df))
        s = map(y -> Symbol(y),y)
        df[!, :date] .= agg.(df[!,:date]);
        df_agg = DataFrames.combine(groupby(df, :date), 
            y .=> fun .=> y);
        return(df_agg)
    end

    """
    finds all ctl files in a directory -> NO subdirectories by default
    findctl(snippet::Union{String,Regex};recurse=false, dir="D:/Wasim/regio/control",suffix=".ctl")
    """
    
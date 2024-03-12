# 
function dfpall(dfs::Vector{DataFrame})
        # df = reduce((left, right) -> 
        # innerjoin(left, right, on = :date,makeunique=true), 
        # dfs)
        df = innerjoin(unique.(dfs, xcol)..., on = xcol, makeunique=true)
        
        y = filter(x->!occursin("date",x), names(df))
        s = map(y -> Symbol(y),y)
        @df df Plots.plot(:date,
                cols(s),
                #yaxis = :log,
                legend = :bottom)
    end
    
    """
    reads, reduces + merges by date
    ds = innerjoin(unique.(dfs, xcol)..., on = xcol, makeunique=true)
    example:
    df = mall(glob(r"qges|qbas|qd"))
    """
    
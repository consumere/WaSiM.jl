# 
function dfp!(regex::AbstractString,dfs::Vector{DataFrame})
        "selects first match and plots..."
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),
            dfs))] |> first
        o = DataFrames.metadata(df)|>collect
        ti = basename(o[1][2])
        if (any(x->occursin("year",x),names(df)))
            s = Symbol.(filter(x->!occursin("year",x),names(df)))
            @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
        else    
        s = Symbol.(filter(x->!occursin("date",x),names(df)))
        @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
        end
    end


    
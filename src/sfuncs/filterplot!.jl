# 
function filterplot!(regex::AbstractString,dfs::Vector{DataFrame})
        "selects first match and add to plot..."
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
        map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
        )] |> first
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        #Symbol(names(df))
        #s = propertynames(df)[Not(end)] #geht auch, aber positionsabhängig
        @df df Plots.plot!(:date,cols(s),legend = :topright)
    end

    # filterplot("win",dfs)
    # filterplot!("qg",dfs)
    # filterplot!("qout",dfs)
    #dfs[5] |>dfp
    #convert(DataFrame,dfs[5]|>DataFrames.metadata|>only)

    
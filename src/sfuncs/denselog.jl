# 
function denselog(regex::AbstractString,dfs::Vector{DataFrame})
        "selects first match and plots..."
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),
            dfs))] |> first
            s = propertynames(df)[Not(end)];
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            @df df density(
                cols(s),
                title=ti,
                yaxis=:log,
                legend = :topright) 
    end

    
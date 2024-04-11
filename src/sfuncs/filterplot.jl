# 
function filterplot(regex::AbstractString,dfs::Vector{DataFrame})
        "selects first match and plots..."
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
        map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
        )] |> first
        dfp(df)
        #indexin(1:length(dfs),
        #map(x->basename(only(DataFrames.metadata(x))[2]),dfs))
    end

    
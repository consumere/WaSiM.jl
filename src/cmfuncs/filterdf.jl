# 
function filterdf(regex::AbstractString,dfs::Vector{DataFrame})
        """
        for namestring, see as dfilter
        selects df from dfvector...
        same as getdf
        """
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),
            dfs))] |> first
        return(df)
    end

    
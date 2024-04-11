# 
function filterdf(str::AbstractString,dfs::Vector{DataFrame})
        """
        for namestring, see as dfilter
        selects df from dfvector...
        same as getdf
        example:
        filterdf("clou",dfs)|>bardfm 
        """
        df = dfs[map(n->occursin(Regex(str,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),
            dfs))] |> first
        
        # filter(n->occursin(Regex(regex,"i"),n),
        # map(x->basename(only(DataFrames.metadata(x))[2]),
        # dfs)
        # )
    end

    
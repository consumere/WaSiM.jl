# 
function getdf(regex::AbstractString,dfs::Vector{Any})
        "selects by regex.."
        #df = getindex(dfs,index)
        
        #y = filter(x->!occursin("date",x),names(df))
        al = filter(x->occursin(Regex(regex,"i"),x),dfs)
        
        df = filter(x->!occursin(r"txt|yrly|nc|png|svg|ftz_0|ftz",x),al)|>
        first |>waread
            
        return(df)
    end
 
    
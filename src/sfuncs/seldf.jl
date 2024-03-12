# 
function seldf(str::String,dfs::Vector{DataFrame})
        str  = join([str,"date"],"|")
        #filter(x->select(occursin(str,x)),dfs)
        return map(k -> k[!,Cols(Regex(str))],dfs)
        
    end

    
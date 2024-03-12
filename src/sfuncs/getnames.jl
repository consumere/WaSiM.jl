# 
function getnames(dfs::DataFrame,rg::String)
        nms = names(dfs)
        m = filter(x->occursin(Regex(rg,"i"),x),nms)
        m = length(m)===1 ? only(Symbol.(m)) : Symbol.(m)
        return(m)
    end

    """
    tst of grouped barplot
    """
    
# 
function regand(v::Vector{Any},xv::Regex)
        """
        here you can put any regex to filter the Vector
        like regand(getnames(dfs),r"tem")
        """
        z = v[(broadcast(x->occursin(xv,x),v))] 
    return(z)
    end


    
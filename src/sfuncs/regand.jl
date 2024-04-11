# 
function regand(v::Vector{Any},xv::Regex)
    z = v[(broadcast(x->occursin(xv,x),v))] 
    return(z)
end

"""
dir_path::String, match::String;file_ending="ctl"
"""

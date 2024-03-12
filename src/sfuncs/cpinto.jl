# 
function cpinto(src::Vector{String}, dst::AbstractString;force=false)
    map(x->cp(x,"$dst/$x";force=force),src)
end


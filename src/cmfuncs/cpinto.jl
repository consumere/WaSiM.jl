# 
function cpinto(src::Vector{String}, dst::AbstractString;force=false)
        """
        mkdir("route-bak")
        cpinto(glob("so_inf"), "route-bak")
        rglob("so_inf")
        force=true will first remove an existing dst.
        """
        map(x->cp(x,"$dst/$x";force=force),src)
    end

    
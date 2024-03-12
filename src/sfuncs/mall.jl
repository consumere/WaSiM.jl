# 
function mall(left::DataFrame, right::DataFrame)
        df = innerjoin(left, right, on = :date,makeunique=true)
        return(df)
    end

    """
    mkdir("route-bak")
    cpinto(glob("so_inf"), "route-bak")
    rglob("so_inf")
    force=true will first remove an existing dst.
    """
    
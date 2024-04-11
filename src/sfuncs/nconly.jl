# 
function nconly()
    z = filter(z -> endswith(z,".nc"),readdir())
    return(z)
end


"""
tff2
"""
